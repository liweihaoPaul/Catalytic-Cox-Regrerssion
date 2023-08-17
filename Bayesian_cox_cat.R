

library("rstan")

#' Set up interval indicators for survival data
#'
#' This function constructs interval indicators for survival data based on the given event times and censoring indicators.
#' It helps in determining which intervals an observation contributes to the risk set and the event set.
#'
#' @param y A vector of observed times (either event or censoring times).
#' @param delta A binary vector indicating event occurrence (1 for event, 0 for censored).
#' @param s A vector of pre-specified cutpoints for interval settings.
#' @param J An integer representing the total number of intervals.
#'
#' @return A list containing:
#'   * `ind.r`: Matrix indicating which intervals an observation is at risk.
#'   * `ind.d`: Matrix indicating which intervals an observation has an event.
#'   * `d`: Vector indicating the total number of events in each interval.
#'   * `ind.r_d`: Matrix indicating risk set minus event set for each observation across intervals.
#'
#' @examples
#' # Provide example data and function calls here.
#'
#' @seealso Related functions or packages, especially the ones used internally.
setting.interval <- function(y, delta, s, J) {
  n <- length(y)
  
  smax	<- max(s)
  
  case0 <- which(delta == 0)
  case1 <- which(delta == 1)
  
  case0yleq <- which(delta == 0 & y <= smax)
  case0ygeq <- which(delta == 0 & y > smax)
  case1yleq <- which(delta == 1 & y <= smax)
  case1ygeq <- which(delta == 1 & y > smax)
  
  
  ind.d <- ind.r <- matrix(0, n, J)
  
  for (i in case1yleq) {
    d.mat.ind	<- min(which(s - y[i] >= 0))
    ind.d[i, d.mat.ind]	<- 1
    ind.r[i, 1:d.mat.ind] <- 1
  }
  
  for (i in case0yleq) {
    cen.j <- min(which(s - y[i] >= 0))
    ind.r[i, 1:cen.j]	<- 1
  }
  
  if (length(union(case1ygeq, case0ygeq)) > 0) {
    ind.r[union(case1ygeq, case0ygeq),]	<- 1
  }
  
  ind.r_d	<- ind.r - ind.d
  
  
  d	<- colSums(ind.d)
  
  list(
    ind.r = ind.r,
    ind.d = ind.d,
    d = d,
    ind.r_d = ind.r_d
  )
}

#' Bayesian Cox Model with Catalytic Priors using Stan
#'
#' This function implements a Bayesian Cox proportional hazards model with catalytic priors using the Stan probabilistic programming language.
#'
#' @param X A matrix of covariates where rows represent observations and columns represent features.
#' @param Y A vector of observed times (either event or censoring times).
#' @param status A binary vector indicating event occurrence (1 for event, 0 for censored).
#' @param tau A hyperparameter for the model. Default is the number of columns in `X`.
#' @param M A scalar defining the number of synthetic samples.
#' @param X.syn A matrix of synthetic covariates.
#' @param T.syn A vector of synthetic event times.
#' @param status.syn A vector indicating systhetic event occurrence (1 for event, 0 for censored).
#' @param h_daga A function representing the hazard rate for synthetic data. By default, it's set to a constant function.
#'
#' @return A `stanfit` object representing the results of the Bayesian model fitting using Stan.
#'
#' @details 
#' This function takes in survival data along with covariates, synthetic covariates, and other parameters to fit a Bayesian Cox proportional hazards model with catalytic priors using Stan. The function returns the results of the Bayesian model fitting.
#'
#' @examples
#' # Provide example data and function calls here.
#'
#' @seealso Related functions or packages, especially the ones used internally. 
#' \code{\link[Stan]{stan_model}}, \code{\link[Stan]{sampling}}
Bayesian_cox_cat_stan <- function(X, Y, status, X.syn, T.syn, status.syn, tau = NULL, h_daga = NULL) {
  # Extract dimensions
  p <- ncol(X)
  n <- nrow(X)
  M <- nrow(X.syn)
  
  # Set default values
  if (is.null(tau)) {
    tau <- p
  }
  
  if (is.null(h_daga)) {
    h_daga <- function(t) {
      sum(status) / mean(Y) / n   
    }
  }
  
  # Hyperparameters for gamma process prior construction 
  eta0 <- 1 
  kappa0 <- 1
  c0 <- 0.2
  s <- c(sort(Y[status == 1]), 2 * max(Y) - max(Y[-which(Y == max(Y))]))
  J <- length(s)
  
  # Calculate intervals
  intv <- setting.interval(Y, status, s, J)
  
  # Compute H.star and alpha0
  H.star <- alpha0 <- numeric(J)
  for (j in 1:J) {
    H.star[j] <- eta0 * s[j] ^ kappa0
    alpha0[j] <- c0 * H.star[j]
  }
  
  hPriorSh <- diff(c(0, alpha0))
  
  # Compute H_daga_Y_star
  H_daga_Y_star <- sapply(1:M, function(index) {
    integrate(Vectorize(h_daga), lower = 0, upper = T.syn[index])$value
  })
  
  # For grouped likelihood
  R_tilde_minus_D_tilde <- intv$ind.r_d
  D_tilde <- intv$ind.d
  
  # Stan model and sampling
  model <- stan_model("bayesian_cox_cat.stan")
  fit <- sampling(
    model, 
    data = list(
      J = J, 
      hPriorSh = hPriorSh, 
      c0 = c0,
      eta0 = eta0, 
      kappa0 = kappa0, 
      n = n,
      p = p, 
      X = X,
      M = M, 
      X_star = X.syn, 
      Y_star = T.syn, 
      status_star = status.syn,
      H_daga_Y_star = H_daga_Y_star, 
      tau_downweight = tau,
      R_tilde_minus_D_tilde = R_tilde_minus_D_tilde, 
      D_tilde = D_tilde
    ),
    chains = 1
  )
  
  return(fit)
}





