# ================================== hanova1 ================================= #
#
#' Posterior sampling for a 1-way hierarchical ANOVA
#'
#' Produces random samples from the posterior distribution of the parameters
#' of a 1-way hierarchical ANOVA model.
#'
#' @param n A numeric scalar.  The size of posterior sample required.
#' @param resp A numeric vector.  Response values.
#' @param fac A vector of class \link{factor} indicating the group from
#'   which the corresponding element of \code{resp} originates.
#'   Must have the same length as \code{resp}.
#' @param prior The log-prior for the parameters of the hyperprior
#'   distribution.  If the user wishes to specify their own prior then
#'   \code{prior} must be an object returned from a call to
#'   \code{\link{set_user_prior}}.
#'   Otherwise, \code{prior} is a character scalar giving the name of the
#'   required in-built prior.
#'   If \code{prior} is not supplied then a default prior is used.
#'   See \strong{Details}.
#' @param hpars A numeric vector.  Used to set parameters (if any) in
#'   an in-built prior.  If \code{prior = cauchy} then \code{hpars} is
#'   a numeric vector of length 2 giving the respective scale parameters
#'   of the half-Cauchy priors for \eqn{\sigma_\alpha} and \eqn{\sigma}.
#' @param param A character scalar.
#'   If \code{param = "trans"} (the default) then the marginal posterior
#'   of hyperparameter vector \eqn{\phi} is reparameterized in terms of
#'   \eqn{log \sigma_\alpha, log \sigma}.
#'   If \code{param = "original"} the original parameterization,
#'   i.e. \eqn{\sigma_\alpha, \sigma} is used.
#'   The former tends to make the optimizations involved in the
#'   ratio-of-uniforms algorithm more stable and to increase the probability
#'   of acceptance, but at the expense of slower function evaluations.
#' @param init A numeric vector. Optional initial estimates sent to
#'   \code{\link[rust]{ru}} in the search for the mode of the posterior
#'   density of (perhaps a subset of) the hyperparameter vector \eqn{\phi}.
#'   If an in-built prior is used then \code{ru} is used to sample from the
#'   marginal posterior density of \eqn{(\sigma_\alpha, \sigma)}, so
#'   \code{init} must have length 2.  Otherwise, \code{init} has length
#'   equal to the argument \code{anova_d} supplied to
#'   \code{\link{set_user_prior}}.
#' @param mu0,sigma0 A numeric scalar.  Mean and standard deviation of a
#'   normal prior for \eqn{\mu}.  Only used if an in-built prior is used
#'   or if \code{anova_d = 2} is supplied in a call to
#'   \code{\link{set_user_prior}} to set a user-defined prior.
#'   The default, \code{sigma0 = Inf}, sets an improper uniform prior
#'   for \eqn{\mu}.
#' @param nrep A numeric scalar.  If \code{nrep} is not \code{NULL} then
#'   \code{nrep} gives the number of replications of the original dataset
#'   simulated from the posterior predictive distribution.
#'   Each replication is based on one of the samples from the posterior
#'   distribution.  Therefore, \code{nrep} must not be greater than \code{n}.
#'   In that event \code{nrep} is set equal to \code{n}.
#' @param ... Optional further arguments to be passed to
#'   \code{\link[rust]{ru}}.
#' @details
#'   Consider \eqn{I} independent experiments in which the \eqn{ni} responses
#'   \strong{\eqn{y}}\eqn{i} from experiment/group \eqn{i} are normally
#'   distributed with mean \eqn{\theta i} and standard deviation \eqn{\sigma}.
#'   The population parameters \eqn{\theta}1, ...,
#'   \eqn{\theta}\eqn{I} are modelled as random samples from a normal
#'   distribution with mean \eqn{\mu} and standard deviation
#'   \eqn{\sigma_\alpha}.  Let \eqn{\phi = (\mu, \sigma_\alpha, \sigma)}.
#'   Conditionally on \eqn{\theta}1, ..., \eqn{\theta}\eqn{I},
#'   \strong{\eqn{y}}1, ..., \strong{\eqn{y}}\eqn{I}
#'   are independent of each other and are independent of \eqn{\phi}.
#'   A \emph{hyperprior} is placed on \eqn{\phi}.
#'   The user can either choose parameter values of a default hyperprior or
#'   specify their own hyperprior using \code{\link{set_user_prior}}.
#'
#'   The \code{\link[rust]{ru}} function in the \code{\link[rust]{rust}}
#'   package is used to draw a random sample from the marginal posterior
#'   of the hyperparameter vector \eqn{\phi}.
#'   Then, conditional on these values, population parameters are sampled
#'   directly from the conditional posterior density of
#'   \eqn{\theta}1, ..., \eqn{\theta}\eqn{I} given \eqn{\phi} and the data.
#'   See the vignette("bang-c-anova-vignette", package = "bang")
#'   for details.
#'
#'   The following priors are specified up to proportionality.
#'
#' \emph{Priors:}
#'
#' \code{prior = "bda"} (the default):
#' \eqn{\pi(\mu, \sigma_\alpha, \sigma) = 1/\sigma,}
#' that is, a uniform prior for \eqn{(\mu, \sigma_\alpha, log \sigma)},
#' for \eqn{\sigma_\alpha > 0} and \eqn{\sigma > 0}.
#' The data must contain at least 3 groups, that is, \code{fac} must have
#' at least 3 levels, for a proper posterior density to be obtained.
#' [See Sections 5.7 and  11.6 of Gelman et al. (2014).]
#'
#' \code{prior = "unif"}:
#' \eqn{\pi(\mu, \sigma_\alpha, \sigma) = 1,}
#' that is, a uniform prior for \eqn{(\mu, \sigma_\alpha, \sigma)},
#' for \eqn{\sigma_\alpha > 0} and \eqn{\sigma > 0}.
#' [See Section 11.6 of Gelman et al. (2014).]
#'
#' \code{prior = "cauchy"}: independent half-Cauchy priors for
#' \eqn{\sigma_\alpha} and \eqn{\sigma} with respective scale parameters
#' \eqn{A_\alpha} and \eqn{A}, that is,
#' \eqn{\pi(\sigma_\alpha, \sigma) =
#'   1 / [(1 + \sigma_\alpha^2 / A_\alpha^2) (1 + \sigma^2 / A^2)].}
#' [See Gelman (2006).]  The scale parameters (\eqn{A_\alpha}, \eqn{A})
#' are specified using \code{hpars} = (\eqn{A_\alpha}, \eqn{A}).
#' The default setting is \code{hpars = c(10, 10).}
#'
#' \emph{Parameterizations for sampling:}
#'
#'  \code{param = "original"} is (\eqn{\mu, \sigma_\alpha, \sigma}),
#'  \code{param = "trans"} (the default) is
#'   \eqn{\phi1 = \mu, \phi2 = log \sigma_\alpha, \phi3 = log \sigma}.
#'
#' @return An object (list) of class \code{"hef"}, which has the same
#'   structure as an object of class "ru" returned from \code{\link[rust]{ru}}.
#'   In particular, the columns of the \code{n}-row matrix \code{sim_vals}
#'   contain the simulated values of \eqn{\phi}.
#'   In addition this list contains the arguments \code{model}, \code{resp},
#'   \code{fac} and \code{prior} detailed above and an \code{n} by \eqn{I}
#'   matrix \code{theta_sim_vals}: column \eqn{i} contains the simulated
#'   values of  \eqn{\theta}\eqn{i}.  Also included are
#'   \code{data = cbind(resp, fac)} and \code{summary_stats} a list
#'   containing: the number of groups \code{I}; the numbers of responses
#'   each group \code{ni}; the total number of observations; the sample mean
#'   response in each group; the sum of squared deviations from the
#'   group means \code{s}; the arguments to \code{hanova1} \code{mu0} and
#'   \code{sigma0}; \code{call}: the matched call to \code{hanova1}.
#'
#' @examples
#' # ======= Late 21st Century Global Temperature Data =======
#'
#' # Extract data for RCP2.6
#' RCP26_2 <- temp2[temp2$RCP == "rcp26", ]
#'
#' # Sample from the posterior under the default `noninformative' flat prior
#' # for (mu, sigma_alpha, log(sigma)).  Ratio-of-uniforms is used to sample
#' # from the marginal posterior for (log(sigma_alpha), log(sigma)).
#' temp_res <- hanova1(resp = RCP26_2[, 1], fac = RCP26_2[, 2])
#'
#' # Plot of sampled values of (sigma_alpha, sigma)
#' plot(temp_res, params = "ru")
#'
#' # Plot of sampled values of (log(sigma_alpha), log(sigma))
#' # (centred at (0,0))
#' plot(temp_res, ru_scale = TRUE)
#'
#' # Plot of sampled values of (mu, sigma_alpha, sigma)
#' plot(temp_res)
#'
#' # Estimated marginal posterior densities of the mean for each GCM
#' plot(temp_res, params = "pop", which_pop = "all", one_plot = TRUE)
#'
#' # Posterior sample quantiles
#' probs <- c(2.5, 25, 50, 75, 97.5) / 100
#' round(t(apply(temp_res$sim_vals, 2, quantile, probs = probs)), 2)
#'
#' # Ratio-of-uniforms information and posterior sample summaries
#' summary(temp_res)
#'
#' # ======= Coagulation time data, from Table 11.2 Gelman et al (2014) =======
#'
#' # With only 4 groups the posterior for sigma_alpha has a heavy right tail if
#' # the default `noninformative' flat prior for (mu, sigma_alpha, log(sigma))
#' # is used.  If we try to sample from the marginal posterior for
#' # (sigma_alpha, sigma) using the default generalized ratio-of-uniforms
#' # runing parameter value r = 1/2 then the acceptance region is not bounded.
#'
#' # Two remedies: reparameterize the posterior and/or increase the value of r.
#'
#' # (log(sigma_alpha), log(sigma)) parameterization, ru parameter r = 1/2
#' coag1 <- hanova1(resp = coagulation[, 1], fac = coagulation[, 2])
#'
#' # (sigma_alpha, sigma) parameterization, ru parameter r = 1
#' coag2 <- hanova1(resp = coagulation[, 1], fac = coagulation[, 2],
#'                param = "original", r = 1)
#'
#' # Values to compare to those in Table 11.3 of Gelman et al (2014)
#' all1 <- cbind(coag1$theta_sim_vals, coag1$sim_vals)
#' all2 <- cbind(coag2$theta_sim_vals, coag2$sim_vals)
#' round(t(apply(all1, 2, quantile, probs = probs)), 1)
#' round(t(apply(all2, 2, quantile, probs = probs)), 1)
#'
#' # Pairwise plots of posterior samples from the group means
#' plot(coag1, which_pop = "all", plot_type = "pairs")
#'
#' # Independent half-Cauchy priors for sigma_alpha and sigma
#' coag3 <- hanova1(resp = coagulation[, 1], fac = coagulation[, 2],
#'                  param = "original", prior = "cauchy", hpars = c(10, 1e6))
#'
#' @seealso The \code{\link[rust]{ru}} function in the \code{\link{rust}}
#'   package for details of the arguments that can be passed to \code{ru} via
#'   \code{hanova1}.
#' @seealso \code{\link{hef}} for hierarchical exponential family models.
#' @seealso \code{\link{set_user_prior}} to set a user-defined prior.
#' @references Gelman, A., Carlin, J. B., Stern, H. S. Dunson, D. B.,
#'  Vehtari, A. and Rubin, D. B. (2014) \emph{Bayesian Data Analysis}.
#'  Chapman & Hall / CRC.
#' @references Gelman, A. (2006) Prior distributions for variance
#'   parameters in hierarchical models. \emph{Bayesian Analysis},
#'   \strong{1}(3), 515-533. \url{https://doi.org/10.1214/06-BA117A}.
#' @export
hanova1 <- function(n = 1000, resp, fac, ..., prior = "default", hpars = NULL,
                    param = c("trans", "original"), init = NULL,
                    mu0 = 0, sigma0 = Inf, nrep = NULL) {
  if (length(resp) != length(fac)) {
    stop("resp and fac must have the same length")
  }
  param <- match.arg(param)
  the_call <- match.call()
  #
  # Create a matrix in which each row contains the response data for
  # a separate level of the factor fac.
  #
  y_mat <- make_resp_matrix(resp, fac)
  #
  # Create list of data and sample summaries
  #
  ds <- hanova1_data(y = y_mat)
  #
  # Create a list that defines the prior and any parameters in the prior
  prior <- check_prior(prior, model = "anova1", hpars, ds$I)
  # Extract anova_d and then delete it from prior.
  anova_d <- prior$anova_d
  prior$anova_d <- NULL
  #
  # If N(mu0, sigma0^2) is being used for mu then add mu0 and sigma0 to
  # the list of values to be passed to the loglikelihood function
  #
  if (anova_d == 2) {
    ds <- c(ds, mu0 = mu0, sigma0 = sigma0)
  }
  # -------------------------- Set up log-posterior -------------------------- #
  #
  if (anova_d == 2) {
    loglik_fn <- two_d_log_marg_lik_anova
  } else {
    loglik_fn <- log_marg_lik_anova
  }
  # Set a model-specific log-posterior
  logpost <- function(x, ds) {
    loglik <- do.call(loglik_fn, c(list(x = x), ds))
    if (is.infinite(loglik)) {
      return(loglik)
    }
    logprior <- do.call(prior$prior, c(list(x), prior[-1]))
    return(loglik + logprior)
  }
  #
  # Extract arguments to be passed to ru() and set default values for
  # trans and rotate if they have not been supplied.
  ru_args <- list(...)
  if (is.null(ru_args$trans)) {
    if (param == "original") {
      ru_args$trans <- "none"
    } else {
      ru_args$trans <- "user"
    }
  }
  if (is.null(ru_args$rotate)) {
    ru_args$rotate <- TRUE
  }
  # Calculate initial estimates
  if (is.null(init)) {
    init <- init1anova(y = y_mat)
    if (anova_d == 2) {
      init <- init[-1]
    }
  } else {
    if (length(init) != 3) {
      stop("init must have length 3")
    }
  }
  #
  # Create list of objects to send to function ru()
  fr_list <- list(trans = ru_args$trans, rotate = ru_args$rotate,
                  param = param, anova_d = anova_d)
  fr <- do.call(one_way_anova_create_ru_list, fr_list)
  #
  # Set ru_args$n_grid and ru_args$ep_bc to NULL just in case they have been
  # specified in ...
  ru_args$n_grid <- NULL
  ru_args$ep_bc <- NULL
  for_ru <- c(list(logf = logpost, ds = ds), fr, list(init = init, n = n),
              ru_args)
  # If we use the transformed parameterization then set the transformation
  # and the log-Jacobian of this transformation.
  if (param == "trans") {
    if (anova_d == 2) {
      phi_to_theta <- one_way_anova_phi_to_theta
    } else {
      phi_to_theta <- three_d_one_way_anova_phi_to_theta
    }
    log_j <- one_way_anova_log_j
    for_ru <- c(for_ru, list(phi_to_theta = phi_to_theta, log_j = log_j))
  }
  #
  res <- do.call(rust::ru, for_ru)
  #
  if (anova_d == 2) {
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
    va <- res$sim_vals[, 1] ^ 2
    ve <- res$sim_vals[, 2] ^ 2
    mv <- t(mapply(mvfun, va, ve))
    mu <- stats::rnorm(n, mv[, 1], sqrt(mv[, 2]))
    res$sim_vals <- cbind(mu, res$sim_vals)
  }
  temp <- one_way_anova_cond_sim(res$sim_vals, ds, n)
  res <- c(res, temp)
  #
  # Add information about the model, the data and the prior
  res$model <- "anova1"
  res$resp <- resp
  res$fac <- fac
  res$data <- cbind(resp, fac)
  res$summary_stats <- ds
  res$prior <- prior
  if (anova_d == 2) {
    res$ru <- 2:3
  } else {
    res$ru <- 1:3
  }
  # Simulate nrep datasets from the posterior predictive distribution
  if (!is.null(nrep)) {
    nrep <- min(nrep, n)
    res$data_rep <- sim_pred_hanova1(res$theta_sim_vals, res$sim_vals, fac,
                                     nrep)
    # Transpose, so that the replicates are in the rows
    res$data_rep <- t(res$data_rep)
  }
  res$call <- the_call
  class(res) <- "hef"
  return(res)
}

# ========================== one_way_anova_cond_sim ========================== #

# Sample from the conditional posterior distribution of the population
# parameters given the hyperparameters and the data

one_way_anova_cond_sim <- function(x, ds, n) {
  mu <- x[, 1]
  va <- x[, 2] ^ 2
  ve <- x[, 3] ^ 2
  theta_sim_vals <- matrix(NA, ncol = ds$I, nrow = n)
  for (i in 1:ds$I) {
    ni <- ds$ni[i]
    div <- va + ve / ni
    mean_alpha <- va * (ds$ybari[i] - mu) / div
    sd_alpha <- sqrt((va * ve / ni) / div)
    theta_sim_vals[, i] <- stats::rnorm(n = n, mean = mean_alpha,
                                        sd = sd_alpha)
  }
  # Add mu to get the group means
  theta_sim_vals <- theta_sim_vals + mu
  colnames(theta_sim_vals) <- paste("theta[",1:ds$I, "]", sep = "")
  return(list(theta_sim_vals = theta_sim_vals))
}

# ============================== hanova1_data ================================ #

hanova1_data <- function(y) {
  #
  # Calculates summary statistics for use in the 1-way random effects ANOVA
  # log-likelihood.
  #
  # Args:
  #   y         : data.  I by max(ni) matrix, where ni is the number of
  #               observations in group i
  #
  # Returns: a list containing the summary statistics
  #
  # Number of groups
  I <- nrow(y)
  # Sample size for each group
  ni <- rowSums(!is.na(y))
  # The sum of the numbers of non-missing observations in the rows of x
  ndot <- sum(ni)
  # Mean of data for each group
  ybari <- rowMeans(y, na.rm=TRUE)
  # Sum of squared deviations from group means
  s <- sum((y - ybari) ^ 2, na.rm = TRUE)
  return(list(I = I, ni = ni, ndot = ndot, ybari = ybari, s = s))
}

# ============================= log_marg_lik_anova =========================== #

log_marg_lik_anova <- function(x, I, ni, ndot, ybari, s) {
  #
  # Calculates the log of the marginal posterior of the parameters
  #   mu : overall mean
  #   sa : standard deviation over groups (sigma_alpha)
  #   se : error standard deviation (sigma)
  # in a 1-way random effects ANOVA model
  #
  # Args:
  #   x        : (mu, sigma_alpha, sigma)
  #   I        : number of groups (levels of the explanatory factor)
  #   ni       : vector of the numbers of observations in group i, i=1,...,I
  #   ndot     : sum of the ni (total number of observations)
  #   ybari    : mean of the observations in group i, i=1,...,I
  #   s        : sum of the squares of the differences between y_ij and ybari
  #
  # Returns:
  #   the value of the log of the marginal posterior
  #
  if (x[2] < 0 | x[3] <= 0) return(-Inf)
  mu <- x[1]
  va <- x[2] ^ 2
  ve <- x[3] ^ 2
  #
  si2 <- va + ve / ni
  log_lik <- ((I - ndot) * log(ve) - s / ve - sum(log(si2)) -
              sum((mu - ybari) ^ 2 / si2)) / 2
  return(log_lik)
}

# ============================= log_marg_lik_anova =========================== #

two_d_log_marg_lik_anova <- function(x, I, ni, ndot, ybari, s, mu0, sigma0) {
  #
  # Calculates the log of the marginal posterior of the parameters
  #   mu : overall mean
  #   sa : standard deviation over groups (sigma_alpha)
  #   se : error standard deviation (sigma)
  # in a 1-way random effects ANOVA model
  #
  # Args:
  #   x        : mu, sigma_alpha, sigma
  #   I        :
  #   ni       :
  #   ndot     :
  #   ybari    :
  #   s        :
  #
  # Returns:
  #   the value of the log of the marginal posterior
  #
  if (x[1] < 0 | x[2] <= 0) return(-Inf)
  va <- x[1] ^ 2
  ve <- x[2] ^ 2
  #
  si2 <- va + ve / ni
  #
  mui_vec <- c(mu0, ybari)
  sigma2i_inv_vec <- 1 / c(sigma0 ^ 2, si2)
  s0 <- sum(sigma2i_inv_vec)
  s1 <- sum(mui_vec * sigma2i_inv_vec)
  s2 <- sum(mui_vec ^ 2 * sigma2i_inv_vec)
  #
  log_lik <- ((I - ndot) * log(ve) - s / ve - sum(log(si2)) -
               log(s0) + s1 ^ 2 / s0 - s2) / 2
  return(log_lik)
}

# ================================= init1anova =============================== #

# Change to REML?  ... to get SEs too? ... for use in trans = "BC" to get lambda.
# Suggests lme4.
init1anova <- function(y) {
  alpha <- rowMeans(y, na.rm = TRUE) # means of each group of data
  mu_hat <- mean(y, na.rm=TRUE)      # overall mean
  # We base the initial estimates of va and ve on the expected values of the
  # mean squares given, for the balanced case, in Section 2. of the
  # Seuss et al document.
  ni <- rowSums(!is.na(y))           # sample size for each group
  I <- nrow(y)                       # number of groups
  MS_E <- sum((y - alpha)^2, na.rm = TRUE) / sum(ni-1)
  ve <- MS_E
  MS_B <- sum(ni * (alpha - mu_hat) ^ 2, na.rm = TRUE) / (I - 1)
  va <- (MS_B - MS_E) / mean(ni)
  init <- c(mu_hat, sqrt(va), sqrt(ve))
  return(init)
}

# ============================== make_resp_matrix ============================ #

make_resp_matrix <- function(resp, fac) {
  #
  # Create a matrix in which each row contains the response data for
  # a separate level of the factor fac.
  #
  # Args:
  #   resp : A numeric vector.  Response values.
  #   fac  : A vector of class factor indicating the group from
  #          which the corresponding element of resp originates.
  #          Must have the same length as resp.
  # Returns:
  #   A numeric matrix with number of rows equal to the number of
  #   levels of the factor fac and number of columns equal to the largest
  #   number of observations across the levels of fac.  Rows are padded with
  #   NAs as necessary.
  #
  fac_names <- unique(fac)
  nlevels <- length(fac_names)
  max_rep <- max(table(fac))
  mat <- matrix(NA, nlevels, max_rep)
  for (i in 1:nlevels) {
    which_rows <- which(fac == fac_names[i])
    mat[i, 1:length(which_rows)] <- resp[which_rows]
  }
  row.names(mat) <- fac_names
  return(mat)
}

# ======================== Create list of arguments for ru() ================= #

one_way_anova_create_ru_list <- function(trans, rotate, param, anova_d) {
  d <- anova_d
  if (anova_d == 2) {
    var_names <- c("sigma[alpha]", "sigma")
  } else {
    var_names <- c("mu", "sigma[alpha]", "sigma")
  }
  if (param == "trans") {
    lower <- rep(-Inf, anova_d)
    upper <- rep(Inf, anova_d)
  } else {
    if (anova_d == 2) {
      lower <- c(0, 0)
      upper <- c(Inf, Inf)
    } else {
      lower <- c(-Inf, 0, 0)
      upper <- c(Inf, Inf, Inf)
    }
  }
  return(list(d = d, lower = lower, upper = upper, var_names = var_names))
}

# Transformation from
# phi = [log(sigma_alpha), log(sigma)]
# theta = (sigma_alpha, sigma)

one_way_anova_phi_to_theta <- function(phi) {
  return(exp(phi))
}

# Log-Jacobian of the transformation from theta to phi, i.e. based on the
# derivatives of phi with respect to theta

one_way_anova_log_j <- function(theta) {
  return(-log(theta[1]) - log(theta[2]))
}

# Transformation from
# phi = [log(sigma_alpha), log(sigma)]
# theta = (sigma_alpha, sigma)

three_d_one_way_anova_phi_to_theta <- function(phi) {
  return(c(phi[1], exp(phi[2:3])))
}

# ============================= sim_pred_hanova1 ============================= #

#' Simulate from a one-way hierarchical ANOVA posterior predictive distribution
#'
#' Simulates \code{nrep} draws from the posterior predictive distribution
#' of the one-way hierarchical ANOVA model described in \code{\link{hanova1}}.
#' This function is called within \code{\link{hanova1}} when the argument
#' \code{nrep} is supplied.
#' @param theta_sim_vals A numeric matrix with \code{length(fac)} columns.
#'   Each row of \code{theta_sim_vals} contains normal means simulated from
#'   their posterior distribution.
#' @param sim_vals A numeric matrix with \code{length(fac)} columns.
#'   Each row of \code{sim_vals} contains normal standard deviations \eqn{\sigma}
#'   simulated from their posterior distribution.
#' @param fac The argument \code{fac} to \code{\link{hanova1}}, that is,
#'   a vector of class \link{factor} indicating group membership.
#' @param nrep A numeric scalar.  The number of replications of the original
#'   dataset simulated from the posterior predictive distribution.
#'   If \code{nrep} is greater than \code{nrow(theta_sim_vals)} then
#'   \code{nrep} is set equal to \code{nrow(theta_sim_vals)}.
#' @return A numeric matrix with \code{nrep} columns.  Each column contains
#'   a draw from the posterior predictive distribution of the number of
#'   successes.
#' @examples
#' RCP26_2 <- temp2[temp2$RCP == "rcp26", ]
#' temp_res <- hanova1(resp = RCP26_2[, 1], fac = RCP26_2[, 2])
#' sim_pred <- sim_pred_hanova1(temp_res$theta_sim_vals, temp_res$sim_vals,
#'                              RCP26_2[, 2], 50)
#' @export
sim_pred_hanova1 <- function(theta_sim_vals, sim_vals, fac, nrep) {
  nrep <- min(nrep, nrow(theta_sim_vals))
  # Extract the first nrep rows from the posterior sample of thetas
  thetas <- theta_sim_vals[1:nrep, , drop = FALSE]
  # Extract the first nrep values in the posterior sample of sigma
  sigma <- sim_vals[1:nrep, "sigma", drop = FALSE]
  # Combine thetas and sigmas
  sigma_and_thetas <- cbind(sigma, thetas)
  # Count the numbers of observations of each level of fac
  n_replicates <- as.numeric(table(fac))
  # Remove any zeros
  n_replicates <- n_replicates[n_replicates > 0]
  # Function to simulate one set of data
  norm_fn <- function(x) {
    mu_vec <- rep(x[-1], times = n_replicates)
    return(stats::rnorm(n = length(mu_vec), mean = mu_vec, sd = x[1]))
  }
  return(apply(sigma_and_thetas, 1, norm_fn))
}
