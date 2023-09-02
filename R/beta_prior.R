# Default priors
# Admissibity of x is checked in the call to the log-likelihood
# so we don't need to check it again here.
# Also useful in case the user supplies their own prior and doesn't check this.

# Prior in Section 5.3 of Gelman et al. (2014) Bayesian Data Analysis.
# Uniform on ( alpha / (alpha + beta), (alpha + beta) ^ (-1/2) )

beta_bda_prior <- function(x) {
  return(-2.5 * log(x[1] + x[2]))
}

# Independent gamma priors on alpha and beta

beta_gamma_prior <- function(x, hpars) {
  return(stats::dgamma(x[1], shape = hpars[1], rate = hpars[2], log = TRUE) +
           stats::dgamma(x[2], shape = hpars[3], rate = hpars[4], log = TRUE))
}

# Default hyperparameters for the exponential prior:
# independent exponentials, each with mean 100.

beta_gamma_hpars <- function() {
  return(c(1, 0.01, 1, 0.01))
}

# Transformation from
# phi   = [log(alpha/beta), log(alpha + beta)] to
# theta = (alpha, beta)

beta_phi_to_theta <- function(phi) {
  apb <- exp(phi[2])
  beta <- apb / (1 + exp(phi[1]))
  return(c(apb - beta, beta))
}

# Log-Jacobian of the transformation from theta to phi, i.e. based on the
# derivatives of phi with respect to theta

beta_log_j <- function(theta) {
  return(-log(theta[1]) - log(theta[2]))
}

# Calculate initial estimates of alpha and beta

beta_init_ests <- function(data, param) {
  # Estimate probabilities
  prob <- data[, 1] / data[, 2]
  # Estimate mean and standard deviation
  # Use method of moments to estimate alpha and beta
  mp <- mean(prob)
  vp <- stats::var(prob)
  if (vp < mp * (1 - mp)) {
    mult <- (mp * (1 - mp) / vp - 1)
    alpha <- mp * mult
    beta <- (1 - mp) * mult
  } else {
    alpha <- 0.1
    beta <- alpha * (1 - mp) / mp
  }
  if (param == "trans") {
    init <- c(log(alpha / beta), log(alpha + beta))
  } else {
    init <- c(alpha, beta)
  }
  return(init)
}

# Create list of arguments for ru()

beta_create_ru_list <- function(param) {
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
