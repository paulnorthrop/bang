context("plot_hef")

my_seed <- 47
my_n <- 5
my_prob <- 0.5

# ------------------------- Beta-binomial model, rat data ------------------- #

rat_res <- hef(model = "beta_binom", data = rat, n = 10)

# No arguments: plot marginal hyperparameter (alpha, beta) posterior
check_NULL <- plot(rat_res, n = my_n, prob = my_prob)
test_that("beta_binom: no arguments gives no error", {
  testthat::expect_identical(check_NULL, NULL)
})

# Parameterization used for sampling from marginal hyperparameter posterior
check_NULL <- plot(rat_res, ru_scale = TRUE, n = my_n, prob = my_prob)
test_that("beta_binom: ru_scale = TRUE givens no error", {
  testthat::expect_identical(check_NULL, NULL)
})

# Same again, but now using params = "ru"
check_NULL <- plot(rat_res, params = "ru", n = my_n, prob = my_prob)
test_that("beta_binom: params = ru", {
  testthat::expect_identical(check_NULL, NULL)
})

# Check that we get an error if we try to plot pairs with only one population
check_error <- try(plot(rat_res, plot_type = "pairs"), silent = TRUE)
test_that("beta_binom: plot_type = pairs gives error if length(which_pop)=1", {
  testthat::expect_identical(class(check_error), "try-error")
})

# If plot_type is supplied ...

# "pairs"
check_NULL <- plot(rat_res, plot_type = "pairs", which_pop = 1:2)
test_that("beta_binom: plot_type = pairs OK when length(which_pop) > 1", {
  testthat::expect_identical(check_NULL, NULL)
})

# "sim"
check_NULL <- plot(rat_res, plot_type = "sim")
test_that("beta_binom: plot_type = sim OK", {
  testthat::expect_identical(check_NULL, NULL)
})

# "dens"
check_NULL <- plot(rat_res, plot_type = "dens")
test_that("beta_binom: plot_type = dens OK", {
  testthat::expect_identical(check_NULL, NULL)
})

# "both"
check_NULL <- plot(rat_res, plot_type = "both")
test_that("beta_binom: plot_type = both OK", {
  testthat::expect_identical(check_NULL, NULL)
})

# "dens" on one plot with legend
check_NULL <- plot(rat_res, plot_type = "dens", which_pop = c(1, 71),
                   one_plot = TRUE)
test_that("beta_binom: plot_type = dens, one plot with legend, OK", {
  testthat::expect_identical(check_NULL, NULL)
})

# "both" on one plot with legend should give an error
check_error <- try(plot(rat_res, plot_type = "both", which_pop = c(1, 71),
                   one_plot = TRUE), silent = TRUE)
test_that("beta_binom: plot_type = both, one plot with legend gives error", {
  testthat::expect_identical(class(check_error), "try-error")
})

# Check that a mis-spelled plot_type gives error
check_error <- try(plot(rat_res, plot_type = "simulated"), silent = TRUE)
test_that("beta_binom: error when plot_type is wrong", {
  testthat::expect_identical(class(check_error), "try-error")
})

# ------------------------- Pump failure data ------------------------------- #

pump_res <- hef(model = "gamma_pois", data = pump, n = my_n)

# "dens"
check_NULL <- plot(pump_res, plot_type = "dens")
test_that("gamma_pois: plot_type = dens OK", {
  testthat::expect_identical(check_NULL, NULL)
})
