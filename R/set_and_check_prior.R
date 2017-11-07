# ========================== set_user_prior ===================================

#' Set a user-defined prior
#'
#' Constructs a user-defined prior distribution for use as the argument
#' \code{prior} in \code{\link{hef}} or \code{\link{hanova1}}.
#'
#' @param prior An R function returning \strong{the log of} the prior density
#'   for of (perhaps a subset of) the hyperparameter vector \eqn{\psi}.
#' @param ... Further arguments giving the names and values of any
#'   parameters involved in the function \code{prior}.
#' @param model A character string.  Abbreviated name of the model:
#'   "beta_binom" for beta-binomial, "gamma_pois" for gamma-Poisson,
#'   "anova1" for 1-way ANOVA.
#' @param anova_d An integer scalar.  Only relevant if \code{model = anova}.
#'   If \code{anova_d = 2} then \code{prior} must return the log-prior
#'   density for the standard deviations \eqn{(\sigma_\alpha, \sigma)}
#'   and a normal prior with mean \code{mu0} and standard deviation
#'   \code{sigma0} is used for \eqn{\mu}.  The values of \code{mu0 = 0} and
#'   \code{sigma0 = Inf} are set in the call to \code{hanova1}, with default
#'   values \code{mu0 = 0} and \code{sigma0 = Inf}.
#'   If \code{anova_d = 3} then \code{prior} must return the log-prior
#'   density for \eqn{(\mu, \sigma_\alpha, \sigma)}.
#' @export
set_user_prior <- function(prior, ..., model = c("beta_binom", "gamma_pois",
                                                 "anova1"), anova_d = 2) {
  if (!is.function(prior)) {
    stop("prior must be a function")
  }
  model <- match.arg(model)
  # Return a list and add additional arguments from ....
  if (model == "anova1") {
    temp <- list(prior = prior, ..., anova_d = anova_d)
  } else {
    temp <- list(prior = prior, ...)
  }
  return(structure(temp, class = "bang_prior", model = model))
}

# ============================ check_prior ====================================

check_prior <- function(prior, model, hpars, n_groups = NULL) {
  # If prior is a character scalar then a default prior is ebing requested
  if (is.character(prior)) {
    prior_name <- prior
    prior <- list()
    # One-way ANOVA
    if (model == "one_way_anova") {
      if (prior_name %in% c("default", "bda") & n_groups < 3) {
        stop("Need >= 3 groups for posterior propriety if bda prior used")
      }
      prior$prior <- switch(prior_name,
                            default = anova1_bda_prior,
                            bda = anova1_bda_prior,
                            unif = anova1_unif_prior,
                            cauchy = anova1_cauchy_prior)
      prior$anova_d <- 2
      if (prior_name == "cauchy") {
        prior$hpars <- hpars
        if (is.null(hpars)) {
          prior$hpars <- anova1_cauchy_hpars()
        }
      }
    }
    # Beta-binomial
    if (model == "beta_binom") {
      prior$prior <- switch(prior_name,
                            default = beta_bda_prior,
                            bda = beta_bda_prior,
                            gamma = beta_gamma_prior)
      if (prior_name == "gamma") {
        prior$hpars <- hpars
        if (is.null(hpars)) {
          prior$hpars <- beta_gamma_hpars()
        }
      }
    }
    # Poisson-gamma
    if (model == "gamma_pois") {
      prior$prior <- switch(prior_name,
                            default = gamma_gamma_prior,
                            gamma = gamma_gamma_prior)
      if (prior_name == "gamma" | prior_name == "default") {
        prior$hpars <- hpars
        if (is.null(hpars)) {
          prior$hpars <- gamma_gamma_hpars()
        }
      }
    }
  } else if (class(prior) == "bang_prior") {
    if (attr(prior, "model") != model) {
      stop("model and model set by set_user_prior() don't match")
    }
  } else {
      stop("A user-defined prior must be set using set_user_prior()")
  }
  #
  return(prior)
}

