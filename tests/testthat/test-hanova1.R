#context("1-way Hierarchical ANOVA log-likelihood")

# Checks that the code to evaluate the log-posterior using the product of
# marginal and conditional posteriors gives the same value (up to an additive
# constant) as a more direct evaluation of the log-posterior.

my_tol <- 1e-5

###############################################################################

# 3D, i.e. marginal posterior for (mu, sigma_alpha, sigma)

three_d_case_fn <- function(x, alpha, resp, fac) {
  #
  # Calculates the value of the one-way ANOVA marginal log-likelihood
  # at x = (mu, sigma_alpha, sigma) based on data matrix data using the
  # obvious coding and non-obvious coding.
  #
  # Args:
  #      x : A numeric vector.  c(mu, sigma_alpha, sigma)
  #  alpha : A numeric vector. c(alpha_1, ..., alpha_I)
  #   resp : A numeric vector.  Response values.
  #    fac : A vector of class \link{factor} indicating the group from
  #          which the corresponding element of \code{resp} originates.
  #          Must have the same length as \code{resp}.
  #
  # Returns:
  #   A list with two components:
  #   obvious : the obvious coding based on lbeta()
  #   clever  : the cleverer, but less obvious coding
  #
  # Data matrix
  y_mat <- make_resp_matrix(resp, fac)
  # Calculate summary statistics
  ds <- hanova1_data(y = y_mat)
  # Function to calculate the marginal log-posterior
  loglik_fn <- log_marg_lik_anova
  # Function for marginal log-posterior for theta = (mu, sigma_alpha, sigma)
  logpost <- function(x, ds) {
    loglik <- do.call(loglik_fn, c(list(x = x), ds))
    return(loglik)
  }
  # Function for conditonal log-posterior of alpha given theta
  cond_logpost <- function(alpha, ds) {
    if (x[2] < 0 | x[3] <= 0) return(-Inf)
    mu <- x[1]
    va <- x[2] ^ 2
    ve <- x[3] ^ 2
    mui <- va * (ds$ybari - mu)
    si1 <- va * ve / ds$ni
    si2 <- va + ve / ds$ni
    cond_mean <- mui / si2
    cond_var <- si1 / si2
    val <- sum(stats::dnorm(alpha - mu, mean = cond_mean, sd = sqrt(cond_var),
                            log = TRUE))
    return(val)
  }
  full_logpost <- function(x, alpha, ds) {
    if (x[2] < 0 | x[3] <= 0) return(-Inf)
    mu <- x[1]
    va <- x[2] ^ 2
    ve <- x[3] ^ 2
    logi <- function(i) {
      y_vec <- na.omit(y_mat[i, ])
      sum(stats::dnorm(y_vec, mean = alpha[i], sd = x[3], log = TRUE))
    }
    term1 <- sum(vapply(1:length(alpha), logi, 0))
    term2 <- sum(stats::dnorm(alpha - mu, mean = 0, sd = x[2], log = TRUE))
    return(term1 + term2)
  }
  res <- list()
  res$indirect <- logpost(x, ds) + cond_logpost(alpha, ds)
  res$direct <- full_logpost(x, alpha, ds)
  return(res$indirect - res$direct)
}

# 100*my_probs% posterior sample quantiles are used as parameter values
my_probs <- c(0.05, 0.5, 0.95)

three_d_test_fn <- function(data, test_string) {
  res <- hanova1(resp = data[, 1], fac = data[, 2])
  theta <- apply(res$sim_vals, 2, quantile, probs = my_probs)
  alpha <- apply(res$theta_sim_vals, 2, quantile, probs = my_probs)
  #
  are_eq <- matrix(NA, 3, 3)
  for (i in 1:3) {
    x <- theta[i, ]
    for (j in 1:3) {
      vec_alpha <- alpha[j, ]
      are_eq[i, j] <- three_d_case_fn(x, vec_alpha, data[, 1], data[, 2])
    }
  }
  diff_are_eq <- diff(range(are_eq))
  testthat::test_that(test_string, {
    testthat::expect_equal(diff_are_eq, 0, tolerance = my_tol)
  })
}

# Test based on the Mid 21st Century CMIP5 global temperature projection data
# for rcp26
my_data <- temp1[temp1$RCP == "rcp26", ]
three_d_test_fn(data = my_data, test_string = "Mid Century CMIP5 Data, rcp26")

# Test based on the Late 21st Century CMIP5 global temperature projection data
# for rcp85
my_data <- temp2[temp2$RCP == "rcp85", ]
three_d_test_fn(data = my_data, test_string = "Late Century CMIP5 Data, rcp85")

# Test based on the coagulation time data
three_d_test_fn(data = coagulation, test_string = "Coagulation Data")

###############################################################################

# 2D, i.e. marginal posterior for (sigma_alpha, sigma)

two_d_case_fn <- function(x, alpha, resp, fac) {
  #
  # Calculates the value of the beta-binomial marginal log-likelihood
  # at x = (alpha, beta) based on data matrix data using the obvious
  # coding and non-obvious coding.
  #
  # Args:
  #      x : A numeric vector.  c(mu, sigma_alpha, sigma)
  #  alpha : A numeric vector. c(alpha_1, ..., alpha_I)
  #   resp : A numeric vector.  Response values.
  #    fac : A vector of class \link{factor} indicating the group from
  #'         which the correspnding element of \code{resp} originates.
  #'         Must have the same length as \code{resp}.
  #
  # Returns:
  #   A list with two components:
  #   obvious : the obvious coding based on lbeta()
  #   clever  : the cleverer, but less obvious coding
  #
  # Data matrix
  y_mat <- make_resp_matrix(resp, fac)
  # Calculate summary statistics
  ds <- hanova1_data(y = y_mat)
  # Add mu0 =0 and sigma0 = Inf so that the prior for mu is uniform
  ds <- c(ds, mu0 = 0, sigma0 = Inf)
  # Function to calculate the marginal log-posterior
  loglik_fn <- two_d_log_marg_lik_anova
  # Function for marginal log-posterior for theta = (mu, sigma_alpha, sigma)
  logpost <- function(x, ds) {
    loglik <- do.call(loglik_fn, c(list(x = x), ds))
    return(loglik)
  }
  # Function for conditonal log-posterior of mu given (sigma_alpha, sigma)
  cond_mu_logpost <- function(x, ds) {
    if (x[2] < 0 | x[3] <= 0) return(-Inf)
    mu <- x[1]
    va <- x[2] ^ 2
    ve <- x[3] ^ 2
    mvfun <- function(va, ve) {
      si2 <- va + ve / ds$ni
      mui_vec <- c(ds$mu0, ds$ybari)
      sigma2i_inv_vec <- 1 / c(ds$sigma0 ^ 2, si2)
      s0 <- sum(sigma2i_inv_vec)
      s1 <- sum(mui_vec * sigma2i_inv_vec)
      m <- s1 / s0
      v <- 1 / s0
      return(c(m, v))
    }
    mv <- mvfun(va, ve)
    val <- sum(stats::dnorm(mu, mean = mv[1], sd = sqrt(mv[2]), log = TRUE))
    return(val)
  }
  # Function for conditonal log-posterior of alpha given theta
  cond_logpost <- function(alpha, ds) {
    if (x[2] < 0 | x[3] <= 0) return(-Inf)
    mu <- x[1]
    va <- x[2] ^ 2
    ve <- x[3] ^ 2
    mui <- va * (ds$ybari - mu)
    si1 <- va * ve / ds$ni
    si2 <- va + ve / ds$ni
    cond_mean <- mui / si2
    cond_var <- si1 / si2
    val <- sum(stats::dnorm(alpha - mu, mean = cond_mean, sd = sqrt(cond_var),
                            log = TRUE))
    return(val)
  }
  full_logpost <- function(x, alpha, ds) {
    if (x[2] < 0 | x[3] <= 0) return(-Inf)
    mu <- x[1]
    va <- x[2] ^ 2
    ve <- x[3] ^ 2
    logi <- function(i) {
      y_vec <- na.omit(y_mat[i, ])
      sum(stats::dnorm(y_vec, mean = alpha[i], sd = x[3], log = TRUE))
    }
    term1 <- sum(vapply(1:length(alpha), logi, 0))
    term2 <- sum(stats::dnorm(alpha - mu, mean = 0, sd = x[2], log = TRUE))
    return(term1 + term2)
  }
  res <- list()
  res$indirect <- logpost(x[-1], ds) + cond_logpost(alpha, ds) +
    cond_mu_logpost(x, ds)
  res$direct <- full_logpost(x, alpha, ds)
  return(res$indirect - res$direct)
}

# 100*my_probs% posterior sample quantiles are used as parameter values
my_probs <- c(0.05, 0.5, 0.95)

two_d_test_fn <- function(data, test_string) {
  res <- hanova1(resp = data[, 1], fac = data[, 2])
  theta <- apply(res$sim_vals, 2, quantile, probs = my_probs)
  alpha <- apply(res$theta_sim_vals, 2, quantile, probs = my_probs)
  #
  are_eq <- matrix(NA, 3, 3)
  for (i in 1:3) {
    x <- theta[i, ]
    for (j in 1:3) {
      vec_alpha <- alpha[j, ]
      are_eq[i, j] <- two_d_case_fn(x, vec_alpha, data[, 1], data[, 2])
    }
  }
  diff_are_eq <- diff(range(are_eq))
  testthat::test_that(test_string, {
    testthat::expect_equal(diff_are_eq, 0, tolerance = my_tol)
  })
}

# Test based on the Mid 21st Century CMIP5 global temperature projection data
# for rcp26
my_data <- temp1[temp1$RCP == "rcp26", ]
two_d_test_fn(data = my_data, test_string = "Mid Century CMIP5 Data, rcp26")

# Test based on the Late 21st Century CMIP5 global temperature projection data
# for rcp85
my_data <- temp2[temp2$RCP == "rcp85", ]
two_d_test_fn(data = my_data, test_string = "Late Century CMIP5 Data, rcp85")

# Test based on the coagulation time data
two_d_test_fn(data = coagulation, test_string = "Coagulation Data")

