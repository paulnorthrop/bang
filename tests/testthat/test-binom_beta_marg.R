#context("Beta-binomial marginal log-likelihood")

# Checks that a quick way to code the beta-binomial marginal log-likelihood,
# that also avoid problems with underflow when alpha and/or beta are very
# large, gives the same value as the obvious way of coding it using the
# lbeta() function.

my_tol <- 1e-5

beta_binom_test_fn <- function(x, data) {
  #
  # Calculates the value of the beta-binomial marginal log-likelihood
  # at x = (alpha, beta) based on data matrix data using the obvious
  # coding and non-obvious coding.
  #
  # Args:
  #      x : numeric vector.  c(alpha, beta)
  #   data : numeric matrix.  Column 1: y, number of successes
  #                           Column 2: n, number of trials
  #
  # Returns:
  #   A list with two components:
  #   obvious : the obvious coding based on lbeta()
  #   clever  : the cleverer, but less obvious coding
  #
  data_list <- binomial_data(data, prior = "dummy")
  y_mat <- data_list$y_mat
  ny_mat <- data_list$ny_mat
  n_mat <- data_list$n_mat
  y <- data[, 1]
  n <- data[, 2]
  res <- list()
  res$obvious <- check_beta_binom_marginal_loglik(x, y, n)
  res$clever <- beta_binom_marginal_loglik(x, y_mat = y_mat, ny_mat = ny_mat,
                                           n_mat = n_mat)
  return(res)
}

# Base the tests on the rat data
# Calculate crude estimates of alpha and beta using the method of moments
x_mom <- beta_init_ests(rat, param = "original")
x <- x_mom

# Test with various powers of the estimates
# Belt and braces!
# At some point (i.e. for large enough alpha and/or beta) the obvious
# code will break down, but it hasn't done this yet!

my_vec <- 2 * (-5:5)
my_vec <- 3 * (-4:4)
for (i in my_vec) {
  for (j in my_vec) {
    x[1] <- x_mom[1] ^ i
    x[2] <- x_mom[2] ^ j
    c_val <- beta_binom_test_fn(x, rat)
    test_string <- paste("alpha =", x[1], "beta =", x[2])
    testthat::test_that(test_string, {
      testthat::expect_equal(c_val$obvious, c_val$clever, tolerance = my_tol)
    })
    x <- rev(x)
    c_val <- beta_binom_test_fn(x, rat)
    test_string <- paste("alpha =", x[1], "beta =", x[2])
    testthat::test_that(test_string, {
      testthat::expect_equal(c_val$obvious, c_val$clever, tolerance = my_tol)
    })
  }
}
