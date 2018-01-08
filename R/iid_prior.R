# Default priors
# Admissibity of x is checked in the call to the log-likelihood
# so we don't need to check it again here.
# Also useful in case the user supplies their own prior and doesn't check this.

# Flat prior

iid_flat_prior <- function(x) {
  return(0)
}
