# For an ANOVA model set the dimension of the prior:
#  1-way: 2 for a prior for (sigma_alpha, sigma) only
#         3 for a prior for (mu, sigma_alpha, sigma)

# ========================== set_user_prior ===================================

#' Set a user-defined prior
#'
#' Describe
#'
#' @param prior An R function.
#' @param model A character string.  Abbreviated name for the
#'   response-mixture combination.
#' @param ... Further arguments giving the names and values of any
#'   parameters involved in the function \code{prior}.
#' @export
set_user_prior <- function(prior, model = c("binom_beta", "pois_gamma",
                                            "anova"), anova_d = 2, ...) {
  if (!is.function(prior)) {
    stop("prior must be a function")
  }
  model <- match.arg(model)
  prior$anova_d <- anova_d
  # Return a list and add additional arguments from ....
  temp <- list(prior = prior, ...)
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
                            cauchy = anova1_norm_cauchy_prior)
      prior$anova_d <- 2
      if (prior_name == "cauchy") {
        prior$hpars <- hpars
        if (is.null(hpars)) {
          prior$hpars <- anova1_cauchy_hpars()
        }
      }
    }
    # Binomial-beta
    if (model == "binom_beta") {
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
    if (model == "pois_gamma") {
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

