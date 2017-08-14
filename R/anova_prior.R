# From page 288 in chapter 11 of Gelman et al. (2014) Bayesian Data Analysis
# (mu, sigma_alpha, log(sigma)) uniform

old_anova1_bda_prior <- function(x) {
  if (x[2] <= 0 | x[3] <= 0) return(-Inf)
  return(-log(x[3]))
}

# (mu, sigma_alpha, sigma) uniform

old_anova1_unif_prior <- function(x) {
  if (x[2] <= 0 | x[3] <= 0) return(-Inf)
  return(0)
}

# Similarly, but with 2D input vector: (sigma_alpha, sigma)

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

old_anova1_norm_cauchy_prior <- function(x, hpars) {
  if (x[2] <= 0 | x[3] <= 0) return(-Inf)
  log_mu <- stats::dnorm(x[1], mean = hpars[1], sd = hpars[2], log = TRUE)
  log_sigma_alpha <- stats::dcauchy(x[2], scale = hpars[3], log = TRUE)
  log_sigma <- stats::dcauchy(x[3], scale = hpars[4], log = TRUE)
  return(log_mu + log_sigma_alpha + log_sigma)
}

# Default hyperparameters for the half-Cauchy prior:

old_anova1_cauchy_hpars <- function() {
  return(c(0, 1e6, 10, 10))
}

# mu normal, sigmas half-Cauchy

anova1_norm_cauchy_prior <- function(x, hpars) {
  if (x[1] <= 0 | x[2] <= 0) return(-Inf)
  log_sigma_alpha <- stats::dcauchy(x[2], scale = hpars[3], log = TRUE)
  log_sigma <- stats::dcauchy(x[3], scale = hpars[4], log = TRUE)
  return(log_sigma_alpha + log_sigma)
}

# Default hyperparameters for the half-Cauchy prior:

anova1_cauchy_hpars <- function() {
  return(c(10, 10))
}
