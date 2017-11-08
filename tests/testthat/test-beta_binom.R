context("hef, beta_binom, in-built prior vs user prior")

my_seed <- 47
my_tol <- 1e-5
my_n <- 10

# ------------------------- Rat tumor data ------------------------------- #

# 1. Default prior

user_prior_fn <- function(x) {
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

test_that("beta-binom: in-built bda = use bda", {
  testthat::expect_equal(rat_res_a$sim_vals, rat_res_b$sim_vals,
                         tolerance = my_tol)
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

test_that("beta-binom: in-built bda = use bda", {
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

test_that("beta-binom: in-built bda = use bda", {
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

test_that("beta-binom: in-built bda = use bda", {
  testthat::expect_equal(rat_res_a$sim_vals, rat_res_b$sim_vals,
                         tolerance = my_tol)
})

