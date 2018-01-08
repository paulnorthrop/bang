# ========================== set_user_prior ===================================

#' Set a user-defined prior
#'
#' Constructs a user-defined prior distribution for use as the argument
#' \code{prior} in \code{\link{hef}} or \code{\link{hanova1}}.
#'
#' @param prior An R function returning \strong{the log of} the prior density
#'   for of (perhaps a subset of) the hyperparameter vector \eqn{\phi}.
#' @param ... Further arguments giving the names and values of any
#'   parameters involved in the function \code{prior}.
#' @param model A character string.  Abbreviated name of the model:
#'   "beta_binom" for beta-binomial and "gamma_pois" for gamma-Poisson
#'   (see \code{\link{hef}}), "anova1" for 1-way ANOVA
#'   (see \code{\link{hanova1}}), "iid" for a random sample from a
#'   univariate distribution (see \code{\link{iid}}).
#' @param anova_d An integer scalar.  Only relevant if \code{model = anova1}.
#'   If \code{anova_d = 2} then \code{prior} must return the log-prior
#'   density for the standard deviations \eqn{(\sigma_\alpha, \sigma)}
#'   and a normal prior with mean \code{mu0} and standard deviation
#'   \code{sigma0} is used for \eqn{\mu}.  The values of \code{mu0 = 0} and
#'   \code{sigma0 = Inf} are set in the call to \code{hanova1}, with default
#'   values \code{mu0 = 0} and \code{sigma0 = Inf}.
#'   If \code{anova_d = 3} then \code{prior} must return the log-prior
#'   density for \eqn{(\mu, \sigma_\alpha, \sigma)}.
#' @param par_names Only relevant if \code{model = "iid"}.
#'   A character vector containing the variable names in the prior, that is,
#'   the names of the parameters about which inference is required.
#'   If these are supplied then they must match exactly (in the same order)
#'   the parameter names of the argument \code{densfun} to \code{\list{iid}}.
#'   Supplying \code{par_names} is optional but it provides an advisible
#'   check that the components of the prior and log-likelihood are matched
#'   correctly.
#' @details For details of the hyperparameters in \eqn{\phi} see the
#'   \strong{Details} section of \code{\link{hef}} for the models
#'   \code{beta_binom} and \code{gamma_pois} and of \code{\link{hanova1}}
#'   for the model \code{anova1}.
#' @return A list of class \code{"bang_prior"}.  Will contain the component
#'   \code{prior}, the user-supplied function to evaluate the log of the prior,
#'   and any arguments supplied in ....
#' @seealso \code{\link{hef}} for hierarchical exponential family models.
#' @seealso \code{\link{hanova1}} for hierarchical one-way analysis of
#'   variance (ANOVA).
#' @examples
#' # User-defined prior, passing parameters
#' # (equivalent to prior = "gamma" with hpars = c(1, 0.01, 1, 0.01))
#' user_prior <- function(x, hpars) {
#'   return(dexp(x[1], hpars[1], log = TRUE) + dexp(x[2], hpars[2], log = TRUE))
#' }
#' user_prior_fn <- set_user_prior(user_prior, hpars = c(0.01, 0.01))
#'
#' #
#' geom_prior <- set_user_prior(dbeta, shape1 = 1, shape2 = 1, log = TRUE,
#'                              model = "iid")
#' @export
set_user_prior <- function(prior, ..., model = c("beta_binom", "gamma_pois",
                                                 "anova1", "iid"),
                           anova_d = 2, par_names = NULL) {
  if (!is.function(prior)) {
    stop("prior must be a function")
  }
  model <- match.arg(model)
  # Return a list and add additional arguments from ....
  if (model == "anova1") {
    temp <- list(prior = prior, ..., anova_d = anova_d)
  } else if (model == "iid") {
    dots_names <- names(list(...))
    prior_args <- formals(prior)
    prior_names <- names(prior_args)
    # Check that all arguments in ... are arguments of prior
    m <- match(dots_names, prior_names)
    if (any(is.na(m))) {
      stop("'...' includes names that are not arguments to 'prior'")
    }
    # Check that all required arguments have a value (apart from the first)
    n_args <- length(prior_args)
    req_names <- lapply(1:n_args, FUN = function(x) is.name(prior_args[[x]]))
    req_names <- prior_names[unlist(req_names)][-1]
    m <- match(req_names, dots_names)
    if (any(is.na(m))) {
      stop("'...' does not include all required arguments to 'prior'")
    }
    if (is.null(par_names)) {
      temp <- list(prior = prior, ...)
    } else {
      temp <- list(prior = prior, ..., par_names = par_names)
    }
  } else {
    temp <- list(prior = prior, ...)
  }
  return(structure(temp, class = "bang_prior", model = model))
}

# ============================ check_prior ====================================

check_prior <- function(prior, model, hpars, n_groups = NULL) {
  # If prior is a character scalar then a default prior is being requested
  if (is.character(prior)) {
    prior_name <- prior
    prior <- list()
    # One-way ANOVA
    if (model == "anova1") {
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
    # A univariate distribution
    if (model == "iid") {
      prior$prior <- switch(prior_name,
                            default = iid_flat_prior,
                            flat = iid_flat_prior)
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

