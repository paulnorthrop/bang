#------------------------------- Binomial-beta -------------------------------#

# Calculate sufficient statistics

poisson_data <- function(data) {
  if (!is.matrix(data) & !is.data.frame(data)) {
    stop("pois_gamma: data must be a matrix or a data frame")
  }
  if (ncol(data) != 2) {
    stop("pois_gamma: data must have 2 columns")
  }
  if (any(data[, 1] < 0)) {
    stop("pois_gamma: the values in column 1 must be non-negative")
  }
  if (any(!is_wholenumber(data[, 1]))) {
    stop("binom_beta: the data must be whole numbers")
  }
  if (any(data[, 2] <= 0)) {
    stop("pois_gamma: the values in column 2 must be positive")
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

pois_gamma_marginal_loglik <- function(x, y, off, sum_y, y_gt_0) {
  if (any(x <= 0)) {
    return(-Inf)
  }
  return(-sum_y * log(x[2]) - sum((x[1] + y) * log(1 + off / x[2])) +
           sum(lgamma(y_gt_0 + x[1]) - lgamma(x[1])))
}

# Sample from the conditional posterior distribution of the population
# parameters given the hyperparameters and the data

pois_gamma_cond_sim <- function(x, data, n_sim) {
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

# Simulate data from the Poisson-gamma model

sim_pois_gamma <- function(n = 1, alpha = 1, beta = 1, off = 1) {
  if (length(off) != 1 & length(off) != n) {
    stop("off must be scalar or a vector of length n")
  }
  theta <- stats::rgamma(n, shape = alpha, rate = beta)
  y <- stats::rpois(n, lambda = off * theta)
  return(y)
}
