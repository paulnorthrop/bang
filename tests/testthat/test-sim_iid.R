context("iid, in-built prior vs user prior")

my_seed <- 47
my_tol <- 1e-5
my_n <- 10

# --------------------------------- Geometric -------------------------------- #

# 1a. User U(0, 1) vs in-built

# Simulate data
x <- rgeom(10, 0.5)

geom_prior <- set_user_prior(dbeta, shape1 = 1, shape2 = 1, log = TRUE,
                             model = "iid", par_names = "prob")
set.seed(my_seed)
res1 <- iid(x, "geometric", prior = geom_prior)
set.seed(my_seed)
res2 <- iid(x, "geometric", prior = "beta")

test_that("iid geometric: in-built user U(0, 1) vs in-built", {
  testthat::expect_equal(res1$sim_vals, res2$sim_vals, tolerance = my_tol)
})

# 1b. User beta(shape1, shape2) vs in-built

# Simulate data
x <- rgeom(10, 0.5)

geom_prior <- set_user_prior(dbeta, shape1 = 2, shape2 = 3, log = TRUE,
                             model = "iid", par_names = "prob")
set.seed(my_seed)
res1 <- iid(x, "geometric", prior = geom_prior)
set.seed(my_seed)
res2 <- iid(x, "geometric", prior = "beta", hpars = c(2, 3))

test_that("iid geometric: in-built user U(0, 1) vs in-built", {
  testthat::expect_equal(res1$sim_vals, res2$sim_vals, tolerance = my_tol)
})
