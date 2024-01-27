

pln.reg <- function(formula, data, n=1500, method = 'BHHH') {
  # This function estimates a Poisson-Lognormal model
  # The formula follows standard R syntax
  # The data input is a dataframe that the model should be estimated on
  # n = the number of Halton draws - Use at least 1000
  # method can take any estimation method allowed in the maxLik library

  mod_df <- model.frame(formula, data)
  X <- model.matrix(formula, data)
  y <- as.numeric(model.response(mod_df))
  x_names <- Names(formula, data)

  p_model <- glm.nb(formula, data)
  start <- as.vector(p_model$coefficients)
  s <- sqrt(log(1/p_model$theta+1))
  start <- append(start, s) # Add starting values for sigma

  poisson_prob <- function(observed, predicted) {
    log_probs <- dpois(observed, predicted, log = TRUE)

    return(log_probs)
  }

  p_poisson_lognormal <- function(p, y, X, n, est_method){

    coefs <- as.array(head(p,-1))
    sigma <- abs(tail(p,1))
    h = halton(n,dim=1)

    probs <- as.vector(rep(0, length(y)))

    i_lognorm <- qnorm(h, mean=0, sd=sigma)

    for (i in i_lognorm){
      coefs_i <- coefs
      coefs_i[1] <- coefs_i[1] + i
      mu <- exp(X %*% coefs_i)
      p_prob <- poisson_prob(y, mu)
      probs <- probs + exp(p_prob)/n
    }

    ll <- sum(log(probs))
    if (est_method == 'bhhh' | est_method == 'BHHH'){
      return(log(probs))
    } else{return(ll)}
  }

  fit <- maxLik(p_poisson_lognormal,
                start = start,
                y = y,
                X = X,
                n = n,
                est_method = method,
                method = method)

  beta_est <- fit$estimate
  beta_pred <- head(beta_est, -1)
  x_names <- append(x_names, 'sigma')
  names(fit$estimate) <- x_names
  fit$estimate['sigma'] <- abs(fit$estimate['sigma'])

  sigma_sq <- fit$estimate['sigma']^2

  fit$formula <- formula
  fit$predictions <- exp(X %*% beta_pred + sigma_sq/2)
  fit$observed <- y
  fit$residuals <- y - fit$predictions
  # Note that this has the predictions, residuals, and observed outcome stored with the model

  return(fit)
}

predict.pln <- function(model, data){
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
  beta_pred <- head(beta_est, -1)

  sigma_sq <- model$estimate['sigma']^2

  predictions <- exp(X %*% beta_pred + sigma_sq/2)
  observed <- y
  residuals <- y - predictions

  pred <- list('prediction'=predictions, 'observed'=observed, 'residuals'=residuals)

  return(pred)
}

