#context("bayesplot")

my_seed <- 47

# ------------------------- Beta-binomial model, rat data ------------------- #

rat_res <- hef(model = "beta_binom", data = rat, n = 10)

# Check that not supplying nrep to hef gives error
check_error <- try(pp_check(rat_res), silent = TRUE)
test_that("pp_check: error when nrep missing", {
  testthat::expect_identical(class(check_error), "try-error")
})

rat_res <- hef(model = "beta_binom", data = rat, n = 10, nrep = 10)

# Check that pp_check returns a ggplot object
check_ggplot <- "ggplot" %in% class(pp_check(rat_res))
test_that("pp_check: check that ggplot object is returned", {
  testthat::expect_identical(check_ggplot, TRUE)
})

# Check that pp_check returns a ggplot object, when raw = TRUE
check_ggplot <- "ggplot" %in% class(pp_check(rat_res, raw = TRUE))
test_that("pp_check: check that ggplot object is returned, when raw = TRUE", {
  testthat::expect_identical(check_ggplot, TRUE)
})
