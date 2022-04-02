#context("hef, gamms_pois, in-built prior vs user prior")

my_seed <- 47
my_tol <- 1e-5
my_n <- 10

# ------------------------- Pump failure data ------------------------------- #

# 1. Default (gamma) prior

user_prior_fn <- function(x, hpars) {
  return(stats::dexp(x[1], hpars[1], log = TRUE) +
           stats::dexp(x[2], hpars[2], log = TRUE))
}
user_prior <- set_user_prior(user_prior_fn, model = "gamma_pois",
                             hpars = c(0.01, 0.01))

# (i) sampling on (rotated) (log(mean), log(alpha + beta)) scale

# In-built
set.seed(my_seed)
pump_res_a <- hef(model = "gamma_pois", data = pump, n = my_n)
# User
set.seed(my_seed)
pump_res_b <- hef(model = "gamma_pois", data = pump, n = my_n,
                 prior = user_prior)

test_that("gamma-pois: in-built gamma = user gamma, param = trans", {
  testthat::expect_equal(pump_res_a$sim_vals, pump_res_b$sim_vals,
                         tolerance = my_tol)
})

# (ii) Default prior, sampling on (alpha, beta) scale

# In-built
set.seed(my_seed)
pump_res_a <- hef(model = "gamma_pois", data = pump, n = my_n, prior = "gamma",
                 param = "original")
# User
set.seed(my_seed)
pump_res_b <- hef(model = "gamma_pois", data = pump, n = my_n,
                 param = "original", prior = user_prior)

test_that("gamma-pois: in-built gamma = user gamma, param = original", {
  testthat::expect_equal(pump_res_a$sim_vals, pump_res_b$sim_vals,
                         tolerance = my_tol)
})

# --------------------------- Simulated data -------------------------------- #

# Simulate one draw from the posterior predictive distribution based on
# the rat data

pump_res <- hef(model = "gamma_pois", data = pump)
sim_data <- sim_pred_gamma_pois(pump_res$theta_sim_vals, pump, 1)
sim_data <- cbind(sim_data, pump[, 2])

# 1. Default (gamma) prior

user_prior_fn <- function(x, hpars) {
  return(stats::dexp(x[1], hpars[1], log = TRUE) +
           stats::dexp(x[2], hpars[2], log = TRUE))
}
user_prior <- set_user_prior(user_prior_fn, model = "gamma_pois",
                             hpars = c(0.01, 0.01))

# (i) sampling on (rotated) (log(mean), log(alpha + beta)) scale

# In-built
set.seed(my_seed)
sim_res_a <- hef(model = "gamma_pois", data = sim_data, n = my_n)
# User
set.seed(my_seed)
sim_res_b <- hef(model = "gamma_pois", data = sim_data, n = my_n,
                  prior = user_prior)

test_that("gamma-pois: sim_data, in-built gamma = user gamma, param = trans", {
  testthat::expect_equal(pump_res_a$sim_vals, pump_res_b$sim_vals,
                         tolerance = my_tol)
})

# (ii) Default prior, sampling on (alpha, beta) scale

# In-built
set.seed(my_seed)
sim_res_a <- hef(model = "gamma_pois", data = sim_data, n = my_n, prior = "gamma",
                  param = "original")
# User
set.seed(my_seed)
sim_res_b <- hef(model = "gamma_pois", data = sim_data, n = my_n,
                  param = "original", prior = user_prior)

test_that("gamma-pois: sim_data, in-built gamma = user gamma, param = original", {
  testthat::expect_equal(pump_res_a$sim_vals, pump_res_b$sim_vals,
                         tolerance = my_tol)
})
