

pll.reg <- function(formula, data, n=1500, method = 'BHHH', max.iters = 200) {
  # This function estimates a Poisson-Lindley-Lognormal mixture model
  # The formula follows standard R syntax
  # The data input is a dataframe that the model should be estimated on
  # n = the number of Halton draws - Use at least 1000
  # method can take any estimation method allowed in the maxLik library
  # max.iters is the maximum number of iterations allowed in the optimization method

  mod_df <- model.frame(formula, data)
  X <- model.matrix(formula, data)
  y <- as.numeric(model.response(mod_df))
  x_names <- Names(formula, data)

  p_model <- glm.nb(formula, data)
  start <- as.vector(p_model$coefficients)
  s <- sqrt(log(1/p_model$theta+1))
  theta <- 0 # starting value for ln(theta)
  start <- append(start, theta) # Add starting values for ln(theta)
  start <- append(start, s) # Add starting values for sigma

  # Define Poisson-Lindley probability Function
  pld.prob <- function(y, mu, theta){
    p <- (theta^2*mu^y*(theta+mu+y+1))/((theta+1)*(mu+theta)^(y+2))
    return(p)
  }

  # define function for adjustment to mean predictions
  PLD.mean.adjustment <- function(theta){
    return((theta+2)/(theta*(theta+1)))
  }

  p_poisson_lindley_lognormal <- function(p, y, X, n, est_method){

    coefs <- as.array(head(p,-2))
    pars <- tail(p,2)
    sigma <- abs(pars[2])
    theta <- exp(pars[1])
    h = halton(n,dim=1)

    probs <- as.vector(rep(0, length(y)))

    i_lognorm <- qnorm(h, mean=0, sd=sigma)

    for (i in i_lognorm){
      coefs_i <- coefs
      coefs_i[1] <- coefs_i[1] + i
      mu <- exp(X %*% coefs_i)
      p_prob <- pld.prob(y, mu, theta)
      probs <- probs + exp(p_prob)/n
    }

    ll <- sum(log(probs))
    if (est_method == 'bhhh' | est_method == 'BHHH'){
      return(log(probs))
    } else{return(ll)}
  }

  fit <- maxLik(p_poisson_lindley_lognormal,
                start = start,
                y = y,
                X = X,
                n = n,
                est_method = method,
                method = method,
                control = list(iterlim = max.iters, printLevel = 2))

  beta_est <- fit$estimate
  beta_pred <- head(beta_est, -2)
  beta.params <- tail(beta_est)
  x_names <- append(x_names, 'ln(theta)')
  x_names <- append(x_names, 'sigma')
  names(fit$estimate) <- x_names
  fit$estimate['sigma'] <- abs(fit$estimate['sigma'])

  sigma_sq <- fit$estimate['sigma']^2
  theta.est <- exp(beta.params[1])

  fit$formula <- formula
  fit$sigma <- fit$estimate['sigma']
  fit$theta <- theta.est
  fit$mean.adjustment <- PLD.mean.adjustment(theta.est)
  fit$predictions <- exp(X %*% beta_pred + sigma_sq/2) * fit$mean.adjustment
  fit$observed <- y
  fit$residuals <- y - fit$predictions
  # Note that this has the predictions, residuals, and observed outcome stored with the model

  return(fit)
}

predict.pll <- function(model, data){
  # This function takes in a Poisson Lognormal model object and dataframe
  # The function returns a dataframe with predictions, observed outcome, and residuals for the data that was input
  poisson_prob <- function(observed, predicted) {
    log_probs <- dpois(observed, predicted, log = TRUE)

    return(log_probs)
  }

  mod_df <- model.frame(model$formula, data)
  X <- model.matrix(model$formula, data)
  y <- as.numeric(model.response(mod_df))

  beta_est <- model$estimate
  beta_pred <- head(beta_est, -2)

  sigma_sq <- model$sigma^2
  adjustment <- model$mean.adjustment

  predictions <- exp(X %*% beta_pred + sigma_sq/2) * adjustment
  observed <- y
  residuals <- y - predictions

  pred <- list('prediction'=predictions, 'observed'=observed, 'residuals'=residuals)

  return(pred)
}

