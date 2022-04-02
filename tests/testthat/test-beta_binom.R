#context("hef, beta_binom, in-built prior vs user prior")

my_seed <- 47
my_tol <- 1e-5
my_n <- 10

# ------------------------- Rat tumor data ------------------------------- #

# 1. Default prior

user_prior_fn <- function(x) {
  if (any(x <= 0)) return(-Inf)
  return(-2.5 * log(x[1] + x[2]))
}
user_prior <- set_user_prior(user_prior_fn)

# (i) sampling on (rotated) (log(mean), log(alpha + beta)) scale

# In-built
set.seed(my_seed)
rat_res_a <- hef(model = "beta_binom", data = rat, n = my_n)
# User
set.seed(my_seed)
rat_res_b <- hef(model = "beta_binom", data = rat, n = my_n,
                 prior = user_prior)

test_that("beta-binom: in-built bda = user bda, param = trans", {
  testthat::expect_equal(rat_res_a$sim_vals, rat_res_b$sim_vals,
                         tolerance = my_tol)
})

# Also check that the summaries are the same

# (i) params = "hyper"

summary_a <- summary(rat_res_a)
summary_b <- summary(rat_res_b)

test_that("beta-binom hyper summary: in-built bda = user bda, param = trans", {
  testthat::expect_equal(summary_a, summary_b, tolerance = my_tol)
})

# (i) params = "pop"

summary_a <- summary(rat_res_a, params = "pop")
summary_b <- summary(rat_res_b, params = "pop")

test_that("beta-binom pop summary: in-built bda = user bda, param = trans", {
  testthat::expect_equal(summary_a, summary_b, tolerance = my_tol)
})

# (ii) Default prior, sampling on (alpha, beta) scale

# In-built
set.seed(my_seed)
rat_res_a <- hef(model = "beta_binom", data = rat, n = my_n,
                 param = "original")
# User
set.seed(my_seed)
rat_res_b <- hef(model = "beta_binom", data = rat, n = my_n,
                 param = "original", prior = user_prior)

test_that("beta-binom: in-built bda = user bda, param = original", {
  testthat::expect_equal(rat_res_a$sim_vals, rat_res_b$sim_vals,
                         tolerance = my_tol)
})

# 2. Gamma prior

user_prior_fn <- function(x, hpars) {
  return(dexp(x[1], hpars[1], log = TRUE) + dexp(x[2], hpars[2], log = TRUE))
}
user_prior <- set_user_prior(user_prior_fn, hpars = c(0.01, 0.01))

# (i) sampling on (rotated) (log(mean), log(alpha + beta)) scale

# In-built
set.seed(my_seed)
rat_res_a <- hef(model = "beta_binom", data = rat, n = my_n, prior = "gamma")
# User
set.seed(my_seed)
rat_res_b <- hef(model = "beta_binom", data = rat, n = my_n,
                 prior = user_prior)

test_that("beta-binom: in-built gamma = user gamma, param = trans", {
  testthat::expect_equal(rat_res_a$sim_vals, rat_res_b$sim_vals,
                         tolerance = my_tol)
})

# (ii) Default prior, sampling on (alpha, beta) scale

# In-built
set.seed(my_seed)
rat_res_a <- hef(model = "beta_binom", data = rat, n = my_n, prior = "gamma",
                 param = "original")
# User
set.seed(my_seed)
rat_res_b <- hef(model = "beta_binom", data = rat, n = my_n,
                 param = "original", prior = user_prior)

test_that("beta-binom: in-built gamma = user gamma, param = original", {
  testthat::expect_equal(rat_res_a$sim_vals, rat_res_b$sim_vals,
                         tolerance = my_tol)
})


# --------------------------- Simulated data -------------------------------- #

# Simulate one draw from the posterior predictive distribution based on
# the rat data

rat_res <- hef(model = "beta_binom", data = rat)
sim_data <- sim_pred_beta_binom(rat_res$theta_sim_vals, rat, 1)
sim_data <- cbind(sim_data, rat[, 2])

# 1. Default prior

user_prior_fn <- function(x) {
  if (any(x <= 0)) return(-Inf)
  return(-2.5 * log(x[1] + x[2]))
}
user_prior <- set_user_prior(user_prior_fn)

# (i) sampling on (rotated) (log(mean), log(alpha + beta)) scale

# In-built
set.seed(my_seed)
sim_res_a <- hef(model = "beta_binom", data = sim_data, n = my_n)
# User
set.seed(my_seed)
sim_res_b <- hef(model = "beta_binom", data = sim_data, n = my_n,
                 prior = user_prior)

test_that("beta-binom: sim_data, in-built bda = user bda, param = trans", {
  testthat::expect_equal(rat_res_a$sim_vals, rat_res_b$sim_vals,
                         tolerance = my_tol)
})

# (ii) Default prior, sampling on (alpha, beta) scale

# In-built
set.seed(my_seed)
sim_res_a <- hef(model = "beta_binom", data = sim_data, n = my_n,
                 param = "original")
# User
set.seed(my_seed)
sim_res_b <- hef(model = "beta_binom", data = sim_data, n = my_n,
                 param = "original", prior = user_prior)

test_that("beta-binom: sim_data, in-built bda = user bda, param = original", {
  testthat::expect_equal(rat_res_a$sim_vals, rat_res_b$sim_vals,
                         tolerance = my_tol)
})

# Check that if init doesn't have length 2 then an error is returned
check_error <- try(hef(model = "beta_binom", data = sim_data, n = my_n,
                       init = 0.1), silent = TRUE)
test_that("beta_binom: error when init has length 1", {
  testthat::expect_identical(class(check_error), "try-error")
})
