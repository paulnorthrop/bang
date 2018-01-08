# Default priors
# Admissibity of x is checked in the call to the log-likelihood
# so we don't need to check it again here.
# Also useful in case the user supplies their own prior and doesn't check this.

### General priors -------------------

# Flat prior

iid_flat_prior <- function(x) {
  return(0)
}

# Beta (shape1, shape2) prior

iid_beta_prior <- function(x, hpars) {
  return(stats::dbeta(x, shape1 = hpars[1], shape2 = hpars[2], log = TRUE))
}

# Default hyperparameters for the beta prior:

iid_beta_hpars <- function() {
  return(c(1, 1))
}

### Distribution-specific priors -------------------

# Geometric

geom_jeffreys_prior <- function(x) {
  return(-log(x) - log(1 - x) / 2)
}
