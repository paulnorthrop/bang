#context("hanova1, in-built prior vs user prior")

my_seed <- 47
my_tol <- 1e-5
my_n <- 10

# ---------------- Late 21st Century Global Temperature Data ---------------- #

# Extract data for RCP2.6
RCP26_2 <- temp2[temp2$RCP == "rcp26", ]

# 1. Default (bda) prior

user_prior_fn <- function(x) {
  if (any(x <= 0)) return(-Inf)
  return(-log(x[2]))
}
user_prior <- set_user_prior(user_prior_fn, model = "anova1")

# (i) sampling on (rotated) (log(mean), log(alpha + beta)) scale

# In-built
set.seed(my_seed)
res26_2_a <- hanova1(resp = RCP26_2[, 1], fac = RCP26_2[, 2], n = my_n)
# User
set.seed(my_seed)
res26_2_b <- hanova1(resp = RCP26_2[, 1], fac = RCP26_2[, 2], n = my_n,
                     prior = user_prior)

test_that("anova1: in-built bda = user bda, param = trans", {
  testthat::expect_equal(res26_2_a$sim_vals, res26_2_b$sim_vals,
                         tolerance = my_tol)
})

# (ii) Default prior, sampling on (alpha, beta) scale

# In-built
set.seed(my_seed)
res26_2_a <- hanova1(resp = RCP26_2[, 1], fac = RCP26_2[, 2], n = my_n,
                     param = "original")
# User
set.seed(my_seed)
res26_2_b <- hanova1(resp = RCP26_2[, 1], fac = RCP26_2[, 2], n = my_n,
                     param = "original", prior = user_prior)

test_that("anova1: in-built bda = user bda, param = original", {
  testthat::expect_equal(res26_2_a$sim_vals, res26_2_b$sim_vals,
                         tolerance = my_tol)
})

# 2. Uniform prior

user_prior_fn <- function(x) {
  if (any(x <= 0)) return(-Inf)
  return(0)
}
user_prior <- set_user_prior(user_prior_fn, model = "anova1")

# (i) sampling on (rotated) (log(mean), log(alpha + beta)) scale

# In-built
set.seed(my_seed)
res26_2_a <- hanova1(resp = RCP26_2[, 1], fac = RCP26_2[, 2], n = my_n,
                     prior = "unif")
# User
set.seed(my_seed)
res26_2_b <- hanova1(resp = RCP26_2[, 1], fac = RCP26_2[, 2], n = my_n,
                     prior = user_prior)

test_that("anova1: in-built unif = user unif, param = trans", {
  testthat::expect_equal(res26_2_a$sim_vals, res26_2_b$sim_vals,
                         tolerance = my_tol)
})

# (ii) Default prior, sampling on (alpha, beta) scale

# In-built
set.seed(my_seed)
res26_2_a <- hanova1(resp = RCP26_2[, 1], fac = RCP26_2[, 2], n = my_n,
                     param = "original", prior = "unif")
# User
set.seed(my_seed)
res26_2_b <- hanova1(resp = RCP26_2[, 1], fac = RCP26_2[, 2], n = my_n,
                     param = "original", prior = user_prior)

test_that("anova1: in-built unif = user unif, param = original", {
  testthat::expect_equal(res26_2_a$sim_vals, res26_2_b$sim_vals,
                         tolerance = my_tol)
})

# 3. Half-Cauchy prior

user_prior_fn <- function(x, hpars) {
  if (any(x <= 0)) return(-Inf)
  log_sigma_alpha <- stats::dcauchy(x[1], scale = hpars[1], log = TRUE)
  log_sigma <- stats::dcauchy(x[2], scale = hpars[2], log = TRUE)
  return(log_sigma_alpha + log_sigma)
}
user_prior <- set_user_prior(user_prior_fn, model = "anova1", hpars = c(10, 20))

# (i) sampling on (rotated) (log(mean), log(alpha + beta)) scale

# In-built
set.seed(my_seed)
res26_2_a <- hanova1(resp = RCP26_2[, 1], fac = RCP26_2[, 2], n = my_n,
                     prior = "cauchy")
# User
set.seed(my_seed)
res26_2_b <- hanova1(resp = RCP26_2[, 1], fac = RCP26_2[, 2], n = my_n,
                     prior = user_prior)

test_that("anova1: in-built unif = user unif, param = trans", {
  testthat::expect_equal(res26_2_a$sim_vals, res26_2_b$sim_vals,
                         tolerance = my_tol)
})

# (ii) Cauchy prior, sampling on (alpha, beta) scale

# In-built
set.seed(my_seed)
res26_2_a <- hanova1(resp = RCP26_2[, 1], fac = RCP26_2[, 2], n = my_n,
                     param = "original", prior = "cauchy", hpars = c(10, 20))
# User
set.seed(my_seed)
res26_2_b <- hanova1(resp = RCP26_2[, 1], fac = RCP26_2[, 2], n = my_n,
                     param = "original", prior = user_prior)

test_that("anova1: in-built unif = user unif, param = original", {
  testthat::expect_equal(res26_2_a$sim_vals, res26_2_b$sim_vals,
                         tolerance = my_tol)
})

# --------------------------- Simulated data -------------------------------- #

# Simulate one draw from the posterior predictive distribution based on
# the rat data

rat_res <- hef(model = "beta_binom", data = rat)
sim_data <- sim_pred_beta_binom(rat_res$theta_sim_vals, rat, 1)
sim_data <- cbind(sim_data, rat[, 2])

RCP26_2 <- temp2[temp2$RCP == "rcp26", ]
res26_2 <- hanova1(resp = RCP26_2[, 1], fac = RCP26_2[, 2])
sim_resp <- sim_pred_hanova1(res26_2$theta_sim_vals, res26_2$sim_vals,
                             RCP26_2[, 2], 1)

# 1. Default (bda) prior

user_prior_fn <- function(x) {
  if (any(x <= 0)) return(-Inf)
  return(-log(x[2]))
}
user_prior <- set_user_prior(user_prior_fn, model = "anova1")

# (i) sampling on (rotated) (log(mean), log(alpha + beta)) scale

# In-built
set.seed(my_seed)
sim_res_a <- hanova1(resp = sim_resp, fac = RCP26_2[, 2], n = my_n)
# User
set.seed(my_seed)
sim_res_b <- hanova1(resp = sim_resp, fac = RCP26_2[, 2], n = my_n,
                     prior = user_prior)

test_that("anova1: sim_data, in-built bda = user bda, param = trans", {
  testthat::expect_equal(sim_res_a$sim_vals, sim_res_b$sim_vals,
                         tolerance = my_tol)
})

# (ii) Default prior, sampling on (alpha, beta) scale

# In-built
set.seed(my_seed)
sim_res_a <- hanova1(resp = sim_resp, fac = RCP26_2[, 2], n = my_n,
                     param = "original")
# User
set.seed(my_seed)
sim_res_b <- hanova1(resp = sim_resp, fac = RCP26_2[, 2], n = my_n,
                     param = "original", prior = user_prior)

test_that("anova1: sim_data, in-built bda = user bda, param = original", {
  testthat::expect_equal(sim_res_a$sim_vals, sim_res_b$sim_vals,
                         tolerance = my_tol)
})


# Check that if init doesn't have length 3 then an error is returned
check_error <- try(hanova1(resp = sim_resp, fac = RCP26_2[, 2], n = my_n,
                           init = 0.1), silent = TRUE)
test_that("anova1: error when init has length 1", {
  testthat::expect_identical(class(check_error), "try-error")
})
