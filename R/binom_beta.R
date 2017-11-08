#------------------------------- Beta-binomial -------------------------------#

# Calculate sufficient statistics

binomial_data <- function(data, prior) {
  if (!is.matrix(data) & !is.data.frame(data)) {
    stop("beta_binom: data must be a matrix or a data frame")
  }
  if (ncol(data) != 2) {
    stop("beta_binom: data must have 2 columns")
  }
  if (any(!is_wholenumber(data))) {
    stop("beta_binom: the data must be whole numbers")
  }
  n <- data[, 2]
  if (any(n <= 0)) {
    stop("beta_binom: the values in data column 2 (n) must be positive")
  }
  y <- data[, 1]
  if (any(y < 0)) {
    stop("beta_binom: the values in data column 1 (y) must be non-negative")
  }
  if (any(y > n)) {
    stop("beta_binom: in data column 1 (y) cannot exceed column 2 (n), rowwise")
  }
  # If a default improper prior is used then check for posterior propriety
  if (is.character(prior)) {
    if (prior == "bda" | prior == "default") {
       if (all(pmin(y, n - y) == 0)) {
         stop("''bda'': improper posterior unless 0 < y < n at least once")
       }
    }
  }
  # Calculate the values needed in beta_binom_marginal_loglik
  ny <- n - y
  tab_y <- table(y[y > 0])
  tab_ny <- table(ny[ny > 0])
  tab_n <- table(n)
  y_vals <- as.numeric(names(tab_y))
  w_y <- as.numeric(tab_y)
  y_mat <- cbind(y_vals, w_y)
  ny_vals <- as.numeric(names(tab_ny))
  w_ny <- as.numeric(tab_ny)
  ny_mat <- cbind(ny_vals, w_ny)
  n_vals <- as.numeric(names(tab_n))
  w_n <- as.numeric(tab_n)
  n_mat <- cbind(n_vals, w_n)
  #
  return(list(y_mat = y_mat, ny_mat = ny_mat, n_mat = n_mat))
}

# Marginal log-likelihood
# posterior for (alpha, beta) not including the prior for (alpha, beta)
# Note: we need to be careful to avoid underflow when either or both
# alpha and beta are very large

beta_binom_marginal_loglik <- function(x, y_mat, ny_mat, n_mat) {
  if (any(x <= 0)) {
    return(-Inf)
  }
  s <- x[1] + x[2]
  mu <- x[1] / s
  f1 <- function(y) {
    return(sum(log(mu + (y - 1:y) / s)))
  }
  t1 <- sum(y_mat[, 2] * vapply(y_mat[, 1], f1, 0))
  f2 <- function(ny) {
    return(sum(log(1 - mu + (ny - 1:ny) / s)))
  }
  t2 <- sum(ny_mat[, 2] * vapply(ny_mat[, 1], f2, 0))
  f3 <- function(n) {
    return(sum(log(1 + (n - 1:n) / s)))
  }
  t3 <- sum(n_mat[, 2] * vapply(n_mat[, 1], f3, 0))
  return(t1 + t2 - t3)
}

# Obvious coding - used only by testthat to test that
# beta_binom_marginal_loglik is correct

check_beta_binom_marginal_loglik <- function(x, y, n) {
  if (any(x <= 0)) {
    return(-Inf)
  }
  return(sum(lbeta(y + x[1], n - y + x[2]) - lbeta(x[1], x[2])))
}

# Sample from the conditional posterior distribution of the population
# parameters given the hyperparameters and the data

beta_binom_cond_sim <- function(x, data, n_sim) {
  alpha <- x[, 1]
  beta <- x[, 2]
  y <- data[, 1]
  n <- data[, 2]
  len_y <- length(y)
  theta_sim_vals <- matrix(NA, ncol = len_y, nrow = n_sim)
  for (i in 1:len_y) {
    theta_sim_vals[, i] <- stats::rbeta(n_sim, alpha + y[i],
                                        beta + n[i] - y[i])
  }
  colnames(theta_sim_vals) <- paste("p[",1:len_y,"]", sep = "")
  return(list(theta_sim_vals = theta_sim_vals))
}

# ----------------------------- sim_beta_binom ------------------------------ #

#' Simulate data from the beta-binomial model
#'
#' Simulates from the beta-binomial model described in \code{\link{hef}}.
#'
#' @param J An integer scalar. The number of groups.
#' @param size A numeric scalar or a numeric vector of length \code{J}.
#'   The number of trials in the groups 1, ..., \code{J}.  If \code{size} is a
#'   scalar then \code{size} is used for all groups.
#' @param alpha,beta Numeric vectors.  The parameters of the
#'   beta(\eqn{\alpha, \beta}) distribution for the binomial success
#'   probabilities \eqn{p1, ..., pJ}.
#' @return A numeric matrix with 2 columns.  The first column, \code{y},
#'   contains the numbers of successes, the second column, \code{n}, the
#'   numbers of trials.
#' @examples
#' # Simulate data that are similar to the rat data
#' sim_data <- sim_beta_binom(J = nrow(rat), size = rat[, "n"], alpha = 2.4,
#'                            beta = 14.3)
#' @export
sim_beta_binom <- function(J = 1, size = 1, alpha = 1, beta = 1) {
  if (length(size) != 1 & length(size) != J) {
    stop("size must be scalar or a vector of length J")
  }
  theta <- stats::rbeta(J, alpha, beta)
  size <- rep_len(size, J)
  y <- mapply(stats::rbinom, size = size, prob = theta, n = 1)
  return(cbind(y, n = size))
}
