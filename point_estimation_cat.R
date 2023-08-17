
# Point estimation --------------------------------------------------------

library(survival)
library(dplyr)

# Risk function -----------------------------------------------------------

#' Identify the Risk Set for a Specific Time in a Survival Dataset
#'
#' This function determines the risk set (observations at risk of experiencing the event) for a specified time point in a survival dataset.
#'
#' @param time_of_interest A scalar representing the time point of interest.
#' @param entry_vector A vector containing entry times for each observation.
#' @param time_vector A vector of event times for each observation.
#' @param event_vector A binary vector indicating event occurrence (1 for event, 0 for censored).
#'
#' @return A vector of indices indicating observations belonging to the risk set at the specified time point.
#'
#' @examples
#' # Provide example data and function calls here.
#'
#' @seealso Related functions or packages, especially the ones used internally.
GetRiskSet <- function(time_of_interest, entry_vector, time_vector, event_vector) {
  return(which((time_of_interest >= entry_vector) & ((time_vector == time_of_interest & event_vector == 1) | (time_vector + 1e-08 > time_of_interest))))
}



# gradient function -------------------------------------------------------
#' Compute the Gradient for Log Partial Likelihood of Cox Proportional Hazards Model
#'
#' This function calculates the gradient vector associated with the log partial likelihood of the Cox proportional hazards model.
#'
#' @param beta A vector representing the coefficients of the model.
#' @param Xs A matrix of covariates where rows represent observations and columns represent features.
#' @param entry A vector of entry times for each observation.
#' @param Ts A vector of event times for each observation.
#' @param event A binary vector indicating event occurrence (1 for event, 0 for censored).
#'
#' @return A vector representing the gradient of the log partial likelihood with respect to beta.
#'
#' @examples
#' # Provide example data and function calls here.
#'
#' @seealso Related functions or packages, especially the ones used internally.
CoxGradient <- function(beta, Xs, entry, Ts, event) {
  p <- ncol(Xs)
  
  gradient <- apply(cbind(Ts, Xs), 1, function(df) {
    df <- matrix(df, nrow = 1)
    ts <- df[, 1]
    xs <- df[, 2:ncol(df)]
    X_risk_set <- Xs[GetRiskSet(ts, entry, Ts, event), ] %>% 
      matrix(ncol = ncol(Xs))
    
    t1 <- t(X_risk_set) %*% exp(X_risk_set %*% beta)
    t2 <- sum(exp(X_risk_set %*% beta))
    
    return(xs - t1 / t2)
  }) %*% event 
  
  return(gradient)
}




#' Compute the Gradient for Log Catalytic Prior for Synthetic Data
#'
#' This function calculates the gradient of the log catalytic prior, specifically for synthetic data without censoring.
#'
#' @param Xsbeta A vector representing the product of the covariates and the coefficients of the model.
#' @param Xs A matrix of covariates where rows represent observations and columns represent features.
#' @param Ts A vector of event times for each observation.
#' @param status.s A vector indicating the status of the synthetic data (typically 1 for event and 0 for censored).
#' @param h_daga A function representing the hazard rate for synthetic data.
#'
#' @return A vector representing the gradient of the log catalytic prior with respect to beta.
#'
#' @examples
#' # Provide example data and function calls here.
#'
#' @seealso Related functions or packages, especially the ones used internally.
CoxstarGradient <- function(Xsbeta, Xs, Ts, status.s, h_daga) {
  M <- nrow(Xs)
  p <- ncol(Xs)
  gd <- rep(0, p)
  
  for (ii in 1:M) {
    gd <- gd + (status.s[ii] - integrate(Vectorize(h_daga), lower = 0, upper = Ts[ii])$value * exp(Xsbeta[ii])) * Xs[ii, ]
  }
  
  return(gd)
}




# Hessian function --------------------------------------------------------
#' Compute the Hessian Matrix for Log Partial Likelihood of Cox Proportional Hazards Model
#'
#' This function calculates the negative Hessian matrix associated with the log partial likelihood of the Cox proportional hazards model.
#'
#' @param beta A vector representing the coefficients of the model.
#' @param Xs A matrix of covariates where rows represent observations and columns represent features.
#' @param entry A vector of entry times for each observation.
#' @param Ts A vector of event times for each observation.
#' @param event A binary vector indicating event occurrence (1 for event, 0 for censored).
#'
#' @return The negative Hessian matrix for the log partial likelihood of the Cox proportional hazards model.
#'
#' @examples
#' # Provide example data and function calls here.
#'
#' @seealso Related functions or packages, especially the ones used internally.
CoxHessian <- function(beta, Xs, entry, Ts, event) {
  p <- ncol(Xs)
  hessian <- matrix(0, p, p)
  uncensored_set <- which(as.numeric(event) == 1)
  
  for (i in uncensored_set) {
    ts <- Ts[i]
    xs <- Xs[i,]
    riskindex <- GetRiskSet(ts, entry, Ts, event)
    
    if (length(riskindex) > 1) {
      X_risk_set <- Xs[riskindex, ]
    } else {
      X_risk_set <- matrix(Xs[riskindex, ], ncol = ncol(Xs))
    }
    
    theta <- exp(X_risk_set %*% beta)
    t1 <- t(X_risk_set) %*% theta
    t2 <- sum(theta)
    
    xxxt <- t(X_risk_set)
    eee <- sapply(1:length(riskindex), function(i) (xxxt[, i] * sqrt(theta[i])))
    
    hessian <- hessian + crossprod(t(eee)) / t2 - t1 %*% t(t1) / t2^2
  }
  
  return(-hessian)  
}



#' Calculate the Hessian Matrix for the Log Catalytic with No Censoring for Synthetic Data
#'
#' This function computes the Hessian matrix for the log catalytic prior, specifically for synthetic data that has no censoring.
#'
#' @param Xsbeta A vector representing the product of the covariates and the coefficients of the model.
#' @param Xs A matrix of covariates (predictor variables) where rows represent observations and columns represent features.
#' @param Ts A vector of event times for each observation.
#' @param h_daga A function representing the hazard rate for synthetic data.
#'
#' @return A Hessian matrix for the log catalytic prior.
#'
#' @examples
#' # Provide example data and function calls here.
#'
#' @seealso Related functions or packages, especially the ones used internally.
CoxstarHessian <- function(Xsbeta, Xs, Ts, h_daga) {
  M <- nrow(Xs)
  p <- ncol(Xs)
  hess <- matrix(0, p, p)
  
  for (ii in 1:M) {
    hess = hess -  integrate(Vectorize(h_daga), lower = 0, upper = Ts[ii])$value* exp(Xsbeta[ii]) * tcrossprod(Xs[ii, ])
  }
  
  return(hess)
}





# RL betahat --------------------------------------------------------------

#' Estimates coefficients using a posterior model with catalytic prior  with Cox's proportional hazards model.
#' Using Newton-Raphson method
#' This function iteratively refines coefficient estimates using gradient and Hessian information.
#' The function leverages both original and synthetic data. Convergence is based on the Euclidean norm of the coefficient difference across iterations.
#'
#' @param initial_values (Optional) Initial coefficient estimates. Default is a vector of zeros.
#' @param TT A vector of original survival times.
#' @param status A vector of original event indicators (1 for event, 0 for censored).
#' @param X A matrix or data frame of original covariates.
#' @param TT.star A vector of synthetic survival times.
#' @param X.star A matrix or data frame of synthetic covariates.
#' @param status.star A vector of synthetic event indicators.
#' @param tau (Optional) A scalar weight. Defaults to the number of columns in `X`.
#' @param h_daga A function representing the h-daga.
#' @param tol Convergence tolerance; the function stops iterating when changes in estimates are below this threshold. Default is 1e-5.
#' @param MAX_ITER Maximum number of iterations. Default is 25.
#'
#' @return A numeric vector of refined coefficient estimates for the Cox proportional hazards model.
#'
#' @examples
#' # Provide example data and function calls here.
#'
#' @seealso Related functions or packages, especially the ones used internally.
betahat_RL_estimation <- function(X,TT, status,
                                   X.star, TT.star, status.star, initial_values = NULL,tau = NULL, 
                                  h_daga=NULL, tol = 1e-5, MAX_ITER = 25) {
  if(is.null(tau)){
    tau<-ncol(X)
  }
  if(is.null(h_daga)){
    h_daga<-function(t){sum(status)/mean(Y)/n}   
  }
  n <- length(TT)
  entry <- rep(0, n)
  if(is.null(initial_values)){
    initial_values=rep(0,ncol(X))
  }
  beta_old <- initial_values
  
  M <- length(TT.star)
  
  Xstarbetaold <- as.numeric(X.star %*% beta_old)
  H_old <- CoxHessian(beta_old, X, entry, TT, status) + tau/M * CoxstarHessian(Xstarbetaold, X.star, TT.star, h_daga)
  gd_old <- CoxGradient(beta_old, X, entry, TT, status) + tau/M * CoxstarGradient(Xstarbetaold, X.star, TT.star,status.star, h_daga)
  beta_new <- beta_old - qr_solve_my(H_old, gd_old) 
  iter <- 1
  diff_new_old <- TRUE
  
  while ((iter < MAX_ITER) & diff_new_old) {
    beta_old <- beta_new
    Xstarbetaold <- as.numeric(X.star %*% beta_old)
    H_old <- CoxHessian(beta_old, X, entry, TT, status) + tau/M * CoxstarHessian(Xstarbetaold, X.star, TT.star, h_daga)
    gd_old <- CoxGradient(beta_old, X, entry, TT, status) + tau/M * CoxstarGradient(Xstarbetaold, X.star, TT.star,status.star, h_daga)
    beta_new <- beta_old - qr_solve_my(H_old, gd_old) 
    iter <- iter + 1
    
    if (is.na(beta_new[1])) {
      warning("NA BLEAK UP")
      return(initial_values + 10000)
    }
    if (norm(as.numeric(beta_new), "2") > 1000000) {
      warning("BLEAK UP")
      return(initial_values + 50000)
    }
    if (norm(as.numeric(beta_old - beta_new), "2") > tol * norm(as.numeric(beta_old), "2")) {
      diff_new_old <- TRUE
    } else {
      diff_new_old <- FALSE
    }
  }
  
  if (iter == MAX_ITER) {
    warning("not converge")
  }
  
  return(as.numeric(beta_new))
}


# Mix betahat -------------------------------------------------------------

#' Estimates the coefficients for a Cox proportional hazards model with mixed data.
#'
#' This function estimates the coefficients for a Cox proportional hazards model by merging and weighting original and synthetic (star) data. If `tau` is not provided, it defaults to a third of the number of columns in `X`.
#'
#' @param X A matrix or data frame of original covariates.
#' @param Y A vector of original survival times.
#' @param status A vector of original event indicators (1 for event, 0 for censored).
#' @param TT.star A vector of synthetic survival times.
#' @param X.star A matrix or data frame of synthetic covariates.
#' @param status.star A vector of synthetic event indicators.
#' @param tau (Optional) A scalar weight for the synthetic data. Defaults to a third of the number of columns in `X`.
#'
#' @return A numeric vector of estimated coefficients for the Cox proportional hazards model using the merged and weighted data.
#' 
#' @examples
#' # Generate some example data
#' X <- matrix(rnorm(100), 50, 2)
#' Y <- rexp(50)
#' status <- sample(0:1, 50, replace = TRUE)
#' TT.star <- rexp(50)
#' X.star <- matrix(rnorm(100), 50, 2)
#' status.star <- sample(0:1, 50, replace = TRUE)
#'
#' # Estimate coefficients
#' betahat_mix_estimation(X, Y, status, TT.star, X.star, status.star)
#' 
#' @seealso \code{\link[survival]{coxph}} for the Cox proportional hazards model function used internally.
betahat_mix_estimation <- function(X, Y, status, X.star,TT.star, status.star, tau = NULL) {
  if (is.null(tau)) {
    tau <- ncol(X) / 3
  }
  
  M <- nrow(X.star)
  X_merge <- rbind(X, X.star)
  Y_merge <- c(Y, TT.star)
  status_merge <- c(status, status.star)
  wt <- c(rep(1, nrow(X)), rep(tau/M, M))
  
  betahat_mix <- as.numeric(coef(coxph(Surv(Y_merge, status_merge) ~ X_merge, weights = wt, ties = "breslow")))
  
  return(betahat_mix)
}




#' Solves a system of linear equations using QR decomposition.
#'
#' This function solves the system of linear equations \(Mx = b\) using the QR decomposition method.
#' If any NA values are encountered in the solution, they are replaced with 0.0001.
#'
#' @param M A matrix of coefficients.
#' @param b A vector representing the constants on the right side of the equations.
#'
#' @return A vector containing the solutions to the system of equations.
#' 
#' 
#' @seealso \code{\link[base]{qr}} for the QR decomposition function used internally.
qr_solve_my <- function(M, b) {
  # Calculate solution using QR decomposition
  s <- qr.coef(qr(M), b)
  
  # Replace any NA values in the solution with 0.0001
  s[is.na(s)] <- 0.0001
  
  return(s)
}

























