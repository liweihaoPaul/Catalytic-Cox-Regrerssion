
# Demo --------------------------------------------------------------------
rm(list=ls())
source("point_estimation_cat.R")
source("Bayesian_cox_cat.R")

# Generation of observed data and synthetic data --------------------------

# Set the random seed for reproducibility
set.seed(123)

# Define the number of features (p) and the sample size (n)
p <- 40
n <- 100

# Define the true beta coefficients
# The first 8 coefficients are specified, while the rest are set to 1 and then divided by sqrt(p)
beta <- c(3, -4, 2, -3, 1, -1, 1, -1, rep(1, p - 8)) / sqrt(p)

# Generate a n x p matrix of random normal values to represent the feature matrix
X <- matrix(rnorm(n * p), n, p)

# Adjust the first two columns of X to be binary (0 or 1) with given probabilities
X[, 1] <- rbinom(n, 1, 0.9)
X[, 2] <- rbinom(n, 1, 0.1)

# Calculate the linear predictor for each observation
Xbeta <- as.numeric(X %*% beta)

# Generate survival times using the exponential distribution with rate 0.5 times exp(Xbeta)
TT <- sapply(1:n, function(i) {
  rexp(1, 0.5 * exp(Xbeta[i]))
})

# Generate censoring times uniformly between 0 and 6
C <- runif(n, 0, 6)

# Determine observed times (Y) and event status
# Y will be the minimum of TT and C
# status will be 1 if the event happened (C > TT), and 0 if the event was censored
Y_status <- t(sapply(1:n, function(aa) {
  c(min(C[aa], TT[aa]), (C[aa] > TT[aa]))
}))

# Separate Y and status into individual vectors
Y <- Y_status[, 1]
status <- Y_status[, 2]


# Generate Synthetic Data for Survival Analysis

# Set the number of synthetic samples
M <- 400

# Function to generate synthetic feature matrix
# Args:
#   X_obs: observed feature matrix
#   M: number of synthetic samples to generate
# Returns:
#   A synthetic feature matrix of size M x p
generate_synthetic_X <- function(X_obs, M) {
  apply(X_obs, 2, function(v) {
    u <- sample(v, size = M, replace = TRUE)
    
    # If the variable is binary, replace some values with random binomial samples
    if (length(unique(v)) <= 2) {
      ind <- sample(c(TRUE, FALSE), M, replace = TRUE)
      u[ind] <- rbinom(sum(ind), 1, 0.5)
    }
    return(u)
  })
}

# Generate synthetic feature matrix using the defined function
X.syn <- generate_synthetic_X(X, M)

# Generate synthetic survival times using exponential distribution
# The rate is determined based on the observed data
T.syn <- rexp(M, rate = sum(status) / mean(Y) / n)

# Set synthetic event status to all 1s (indicating all events occurred)
status.syn <- rep(1, M)



# perform point estimation using catalytic prior with default tau --------------------------------

betahat_RL_estimation(X, Y, status,  X.syn, T.syn,status.syn)
betahat_mix_estimation(X, Y, status,  X.syn, T.syn,status.syn)


# perform Bayesian sampling with stan with default setting -------------------------------------

Bayesian_fit<- Bayesian_cox_cat_stan(X, Y, status, X.syn, T.syn, status.syn)
summary(Bayesian_fit)$summary[1:p,]


