# Default priors
# Admissibity of x is checked in the call to the log-likelihood
# so we don't need to check it again here.
# Also useful in case the user supplies their own prior and doesn't check this.

# Independent gamma priors on alpha and beta

gamma_gamma_prior <- function(x, hpars) {
  return(stats::dgamma(x[1], shape = hpars[1], rate = hpars[2], log = TRUE) +
           stats::dgamma(x[2], shape = hpars[3], rate = hpars[4], log = TRUE))
}

# Default hyperparameters for the exponential prior:
# independent exponentials, each with mean 100.

gamma_gamma_hpars <- function() {
  return(c(1, 0.01, 1, 0.01))
}

# Transformation from
# phi   = [log(alpha/beta), log(beta)] to
# theta = (alpha, beta)

gamma_phi_to_theta <- function(phi) {
  beta <- exp(phi[2])
  return(c(beta * exp(phi[1]), beta))
}

# Log-Jacobian of the transformation from theta to phi, i.e. based on the
# derivatives of phi with respect to theta

gamma_log_j <- function(theta) {
  return(-log(theta[1]) - log(theta[2]))
}

# Calculate initial estimates of alpha and beta

gamma_init_ests <- function(data, param) {
  # Estimate mean and standard deviation
  # Use method of moments to estimate alpha and beta
  rate <- data[, 1] / data[, 2]
  m <- mean(rate)
  v <- stats::var(rate)
  beta <- m / v
  alpha <- beta * m
  if (param == "trans") {
    init <- c(log(alpha / beta), log(beta))
  } else {
    init <- c(alpha, beta)
  }
  return(init)
}

# Create list of arguments for ru()

gamma_create_ru_list <- function(param) {
  d <- 2L
  var_names <- c("alpha", "beta")
  if (param == "trans") {
    lower <- c(-Inf, -Inf)
    upper <- c(Inf, Inf)
  } else {
    lower <- c(0, 0)
    upper <- c(Inf, Inf)
  }
  return(list(d = d, lower = lower, upper = upper, var_names = var_names))
}
