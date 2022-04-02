#context("set_user_prior")

# Non-function prior
check_error <- try(set_user_prior("bda"), silent = TRUE)
test_that("set_user_prior: character prior gives error", {
  testthat::expect_identical(class(check_error), "try-error")
})

