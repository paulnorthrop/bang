#------------------------------- gamma-Poisson -------------------------------#

# Calculate sufficient statistics

poisson_data <- function(data) {
  if (!is.matrix(data) & !is.data.frame(data)) {
    stop("gamma_pois: data must be a matrix or a data frame")
  }
  if (ncol(data) != 2) {
    stop("gamma_pois: data must have 2 columns")
  }
  if (any(data[, 1] < 0)) {
    stop("gamma_pois: the values in column 1 must be non-negative")
  }
  if (any(!is_wholenumber(data[, 1]))) {
    stop("gamma_pois: the data must be whole numbers")
  }
  if (any(data[, 2] <= 0)) {
    stop("gamma_pois: the values in column 2 must be positive")
  }
  y <- data[, 1]
  y_gt_0 <- y[y > 0]
  sum_y <- sum(y_gt_0)
  return(list(y = data[, 1], off = data[, 2], sum_y = sum_y, y_gt_0 = y_gt_0))
}

# Marginal log-likelihood
# posterior for (alpha, beta) not including the prior for (alpha, beta)
# Note: we need to be careful to avoid underflow when either or both
# alpha and beta are very large

gamma_pois_marginal_loglik <- function(x, y, off, sum_y, y_gt_0) {
  if (any(x <= 0)) {
    return(-Inf)
  }
  return(-sum_y * log(x[2]) - sum((x[1] + y) * log(1 + off / x[2])) +
           sum(lgamma(y_gt_0 + x[1]) - lgamma(x[1])))
}

# Sample from the conditional posterior distribution of the population
# parameters given the hyperparameters and the data

gamma_pois_cond_sim <- function(x, data, n_sim) {
  alpha <- x[, 1]
  beta <- x[, 2]
  y <- data[, 1]
  off <- data[, 2]
  len_y <- length(y)
  theta_sim_vals <- matrix(NA, ncol = len_y, nrow = n_sim)
  for (i in 1:len_y) {
    theta_sim_vals[, i] <- stats::rgamma(n_sim, shape = alpha + y[i],
                                         rate = beta + off[i])
  }
  colnames(theta_sim_vals) <- paste("lambda[",1:len_y,"]", sep = "")
  return(list(theta_sim_vals = theta_sim_vals))
}

# --------------------------- sim_pred_gamma_pois --------------------------- #

#' Simulate from a gamma-Poisson posterior predictive distribution
#'
#' Simulates \code{nrep} draws from the posterior predictive distribution
#' of the beta-binomial model described in \code{\link{hef}}.
#' This function is called within \code{\link{hef}} when the argument
#' \code{nrep} is supplied.
#' @param theta_sim_vals A numeric matrix with \code{nrow(data)} columns.
#'   Each row of \code{theta_sim_vals} contains binomial success probabilities
#'   simulated from their posterior distribution.
#' @param data A 2-column numeric matrix: the numbers of successes in column 1
#'   and the corresponding numbers of trials in column 2.
#' @param nrep A numeric scalar.  The number of replications of the original
#'   dataset simulated from the posterior predictive distribution.
#'   If \code{nrep} is greater than \code{nrow(theta_sim_vals)} then
#'   \code{nrep} is set equal to \code{nrow(theta_sim_vals)}.
#' @return A numeric matrix with \code{nrep} columns.  Each column contains
#'   a draw from the posterior predictive distribution of the number of
#'   successes.
#' @examples
#' pump_res <- hef(model = "gamma_pois", data = pump)
#' pump_sim_pred <- sim_pred_gamma_pois(pump_res$theta_sim_vals, pump, 50)
#' @export
sim_pred_gamma_pois <- function(theta_sim_vals, data, nrep) {
  nrep <- min(nrep, nrow(theta_sim_vals))
  # Extract the first nrep rows from the posterior sample of thetas
  thetas <- theta_sim_vals[1:nrep, , drop = FALSE]
  # Function to simulate one set of data
  pois_fn <- function(x) {
    return(stats::rpois(n = length(data[, 2]), lambda = x * data[, 2]))
  }
  return(apply(thetas, 1, pois_fn))
}
