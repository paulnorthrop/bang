# From page 288 in chapter 11 of Gelman et al. (2014) Bayesian Data Analysis
# (mu, sigma_alpha, log(sigma)) uniform

anova1_bda_prior <- function(x) {
  if (x[1] <= 0 | x[2] <= 0) return(-Inf)
  return(-log(x[2]))
}

# (mu, sigma_alpha, sigma) uniform

anova1_unif_prior <- function(x) {
  if (x[1] <= 0 | x[2] <= 0) return(-Inf)
  return(0)
}

# mu normal, sigmas half-Cauchy

anova1_cauchy_prior <- function(x, hpars) {
  if (x[1] <= 0 | x[2] <= 0) return(-Inf)
  log_sigma_alpha <- stats::dcauchy(x[1], scale = hpars[1], log = TRUE)
  log_sigma <- stats::dcauchy(x[2], scale = hpars[2], log = TRUE)
  return(log_sigma_alpha + log_sigma)
}

# Default hyperparameters for the half-Cauchy prior:

anova1_cauchy_hpars <- function() {
  return(c(10, 10))
}
