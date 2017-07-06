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
set_user_prior <- function(prior, model = c("binom_beta", "pois_gamma"), ...) {
  if (!is.function(prior)) {
    stop("prior must be a function")
  }
  model <- match.arg(model)
  # Return a list and add additional arguments from ....
  temp <- list(prior = prior, ...)
  return(structure(temp, class = "bang_prior", model = model))
}

# ============================ check_prior ====================================

check_prior <- function(prior, model, hpars) {
  if (is.character(prior)) {
    prior_name <- prior
    prior <- list()
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
    if (model == "pois_gamma") {
      prior$prior <- switch(prior_name,
                            default = gamma_gamma_prior,
                            exp = gamma_gamma_prior)
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
  return(prior)
}

