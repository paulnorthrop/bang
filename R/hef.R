#' Hierarchical Exponential Family Model
#'
#' Samples from the posterior distribution of the parameters of
#' certain hierarchical exponential family models.  Conditional on
#' population-specific parameter vectors \eqn{\theta}1, ..., \eqn{\theta}J
#' the observed \emph{response} data \eqn{y}1, ..., \eqn{y}J within each
#' population are modelled as random samples from a distribution in an
#' exponential family. The population parameters \eqn{\theta}1, ...,
#' \eqn{\theta}J are modelled as random samples from a common
#' \emph{population distribution}, chosen to be conditionally conjugate
#' to the response distribution, with \emph{hyperparameter} vector \eqn{\phi}.
#' Conditionally on \eqn{\theta}1, ..., \eqn{\theta}J, \eqn{y}1, ..., \eqn{y}J
#' are independent of each other and are independent of \eqn{\psi}.
#' A \emph{hyperprior} is placed on \eqn{\psi}.  The user can either
#' choose parameter values of a default hyperprior or specify their own
#' hyperprior.
#'
#' @param n An integer scalar.  The size of the posterior sample required.
#' @param model A character string.  Abbreviated name for the
#'   response-population distribution combination.
#' @param data A numeric matrix.  The format depends on \code{model}.
#'   See \strong{Details}.
#' @param prior The log-prior for the parameters of the hyperprior
#'   distribution.  If the user wishes to specify their own prior then
#'   \code{prior} must be an object returned from a call to
#'   \code{\link{set_user_prior}}.
#'   If \code{prior} is not supplied then a default prior is used.
#'   See \strong{Details}.
#' @param hpars A numeric vector.  Used to set parameters (if any) in
#'   a default prior.
#' @param param A character scalar.
#'   If \code{param = "trans"} (the default) then the marginal posterior
#'   for of hyperparameter vector \eqn{psi} is parameterised in a way
#'   designed to improve the efficiency of sampling from this posterior.
#'   If \code{param = "original"} the original parameterisation is used.
#'   The former tends to make the optimisations involved in the
#'   ratio-of-uniforms algorithm more stable and to increase the probability
#'   of acceptance, but at the expense of slower function evaluations.
#' @param ... Optional further arguments to be passed to
#'   \code{\link[rust]{ru}}.
#' @details
#'   The \code{\link[rust]{ru}} function is used to draw a random sample
#'   from the marginal posterior of the hyperparameter vector \eqn{\psi}.
#'   Then, conditional on these values, population parameters are sampled
#'   directly from the conditional posterior density of
#'   \eqn{\theta}1, ..., \eqn{\theta}J given \eqn{\psi} and the data.
#'
#'   We outline each \code{model}, specify the format of the
#'   \code{data}, give the default (log-)priors (up to an additive constant)
#'   and detail the choices of ratio-of-uniforms parameterisation
#'   \code{param}.
#'
#' \strong{Binomial-beta:} For \eqn{j = 1, ..., J},
#'   \eqn{Yj | \thetaj} are i.i.d binomial\eqn{(nj, \thetaj)},
#'   where \eqn{\thetaj} is the probability of success in group \eqn{j}
#'   and \eqn{nj} is the number of trials in group \eqn{j}.
#'   \eqn{\thetaj} are i.i.d. beta\eqn{(\alpha, \beta)}, so
#'   and \eqn{\psi = (\alpha, \beta)}.
#'   \code{data} is a 2-column matrix: the numbers of successes in column 1
#'   and the corresponding numbers of trials in column 2.
#'
#' \emph{Priors:}
#'
#' \code{prior = "bda"} (the default):
#' \eqn{log \pi(\alpha, \beta) = - 2.5 log(\alpha + \beta),
#'   \alpha > 0, \beta > 0.} [See Section 5.3 of Gelman et al. (2013).]
#'
#' \code{prior = "gamma"}: independent gamma priors on \eqn{\alpha}
#' and \eqn{\beta}, i.e.
#' \eqn{log \pi(\alpha, \beta) =
#'   (s1 - 1)log\alpha - r1 \alpha +
#'   (s2 - 1)log\beta - r2 \beta,  \alpha > 0, \beta > 0.}
#' where the respective shape (\eqn{s1}, \eqn{s2}) and rate
#' (\eqn{r1}, \eqn{r2}) parameters are specified using
#' \code{hpars} = \eqn{(s1, r1, s2, r2)}.
#'
#' \emph{Parameterisations for sampling:}
#'
#'  \code{param = "original"} is (\eqn{\alpha, \beta}),
#'  \code{param = "trans"} (the default) is
#' \eqn{\phi1 = logit(\alpha/(\alpha+\beta)) = log(\alpha/\beta),
#'   \phi2 = log(\alpha+\beta)}.
#' See Section 5.3 of Gelman et al. (2013).
#'
#' \strong{Poisson-gamma:} For \eqn{j = 1, ..., J},
#'   \eqn{Yj | \thetaj} are i.i.d Poisson\eqn{(ej\thetaj)},
#'   where
#'   \eqn{ej} is the \emph{exposure} in group \eqn{j}, based on the
#'   total length of observation time and/or size of the population at
#'   risk of the event of interest and \eqn{\thetaj} is the mean number
#'   of events per unit of exposure.
#'   \eqn{\thetaj} are i.i.d. gamma\eqn{(\alpha, \beta)}, so
#'   \eqn{\psi = (\alpha, \beta)}.
#'   \code{data} is a 2-column matrix: the counts \eqn{yj} in column 1
#'   and the corresponding exposures \eqn{ej} in column 2.
#'
#' \emph{Priors:}
#'
#' \code{prior = "gamma"} (the default): independent gamma priors
#'   on \eqn{\alpha} and \eqn{\beta}, i.e.
#' \eqn{log \pi(\alpha, \beta) =
#'   (s1 - 1)log\alpha - r1 \alpha +
#'   (s2 - 1)log\beta - r2 \beta,  \alpha > 0, \beta > 0.}
#' where the respective shape (\eqn{s1}, \eqn{s2}) and rate
#' (\eqn{r1}, \eqn{r2}) parameters are specified using
#' \code{hpars} = \eqn{(s1, r1, s2, r2)}.
#'
#' \emph{Parameterisations for sampling:}
#'
#'  \code{param = "original"} is (\eqn{\alpha, \beta}),
#'  \code{param = "trans"} (the default) is
#' \eqn{\phi1 = log(\alpha/\beta), \phi2 = log(\beta).}
#' @return An object (list) of class \code{"hef"}, which has the same
#'   structure as an object of class "ru" returned from \code{\link[rust]{ru}}.
#'   In particular, the columns of the \code{n}-row matrix \code{sim_vals}
#'   contain the simulated vaues of \eqn{psi}.
#'   In addition this list contains the arguments \code{model}, \code{data},
#'   \code{prior} detailed above and an \code{n} by \eqn{J} matrix
#'   \code{theta_sim_vals}: column j contains the simulated values of
#'   (the first component of) \eqn{\thetaj}.  If \eqn{\thetaj} has a
#'   second component then the simulated values are given in
#'   \code{theta_2_sim_vals}.
#' @references Gelman, A., Carlin, J. B., Stern, H. S. Dunson, D. B.,
#'  Vehtari, A. and Rubin, D. B. (2013) \emph{Bayesian Data Analysis}.
#'  Chapman & Hall / CRC.
#'   \url{http://www.stat.columbia.edu/~gelman/book}

#' @examples
#' ############################ Binomial-Beta #################################
#'
#' # ------------------------- Rat tumor data ------------------------------- #
#'
#' # Default prior, sampling on (rotated) (log(mean), log(alpha + beta)) scale
#' rat_res <- hef(model = "binom_beta", data = rat)
#' plot(rat_res)
#' plot(rat_res, ru_scale = TRUE)
#' summary(rat_res)
#'
#' # Default prior, sampling on (rotated) (alpha, beta) scale
#' rat_res <- hef(model = "binom_beta", data = rat, param = "original")
#' plot(rat_res)
#' plot(rat_res, ru_scale = TRUE)
#' summary(rat_res)
#'
#' # To produce a plot akin to Figure 5.3 of Gelman et al. (2013) we
#' # (a) Use the same prior for (alpha, beta)
#' # (b) Don't use axis rotation (rotate = FALSE)
#' # (c) Plot on the scale used for ratio-of-uniforms sampling (ru_scale = TRUE)
#' # (d) Note that the mode is relocated to (0, 0) in the plot
#' rat_res <- hef(model = "binom_beta", data = rat, rotate = FALSE)
#' plot(rat_res, ru_scale = TRUE)
#' # This is the estimated location of the posterior mode
#' rat_res$f_mode
#'
#' # User prior, with parameters (equivalent to prior = "exp")
#' user_prior <- function(x, hpars) {
#'   return(dexp(x[1], hpars[1], log = TRUE) + dexp(x[2], hpars[2], log = TRUE))
#' }
#' user_prior <- set_user_prior(user_prior, hpars = c(0.01, 0.01))
#' rat_res <- hef(model = "binom_beta", data = rat, prior = user_prior)
#' plot(rat_res)
#' summary(rat_res)
#'
#' ############################ Poisson-gamma #################################
#'
#' # ------------------------ Pump failure data ------------------------------ #
#'
#' pump_res <- hef(model = "pois_gamma", data = pump)
#' plot(pump_res)
#' plot(pump_res, ru_scale = TRUE)
#' summary(pump_res)
#' @export
hef <- function(n = 1000, model = c("binom_beta", "pois_gamma"),
                data, prior = "default",
                hpars = NULL, param = c("trans", "original"), ...) {
  model <- match.arg(model)
  param <- match.arg(param)
  #
  # Calculate the data summaries required in the posterior
  # and check for posterior propriety, where possible
  ds <- switch(model,
               binom_beta = binomial_data(data, prior),
               pois_gamma = poisson_data(data))
  #
  # Create a list that defines the prior and any parameters in the prior
  prior <- check_prior(prior, model, hpars)
  #
  # Set the function to calculate the log-likelihood
  loglik_fn <- switch(model,
                      binom_beta = binom_beta_marginal_loglik,
                      pois_gamma = pois_gamma_marginal_loglik)
  #
#  gs_init <- gs_alpha(data)
  #
  # Set a model-specific log-posterior
  logpost <- function(x, ds) {
    loglik <- do.call(loglik_fn, c(list(x = x), ds))
    if (is.infinite(loglik)) {
      return(loglik)
    }
    logprior <- do.call(prior$prior, c(list(x), prior[-1]))
    return(loglik + logprior)
  }
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
  init <- switch(model,
                 binom_beta = beta_init_ests(data, param = param),
                 pois_gamma = gamma_init_ests(data, param = param))
  #
  # Create list of objects to send to function ru()
  fr_list <- list(model = model, trans = ru_args$trans,
                  rotate = ru_args$rotate, param = param)
  fr <- switch(model,
               binom_beta = do.call(beta_create_ru_list, fr_list),
               pois_gamma = do.call(gamma_create_ru_list, fr_list))
  #
  # Set ru_args$n_grid and ru_args$ep_bc to NULL just in case they have been
  # specified in ...
  ru_args$n_grid <- NULL
  ru_args$ep_bc <- NULL
  for_ru <- c(list(logf = logpost, ds = ds), fr, list(init = init, n = n),
              ru_args)
  #
  # If we use the transformed parameterisation then set the transformation
  # and the log-Jacobian of this transformation.
  if (param == "trans") {
    phi_to_theta <- switch(model,
                           binom_beta = beta_phi_to_theta,
                           pois_gamma = gamma_phi_to_theta)
    log_j <- switch(model,
                    binom_beta = beta_log_j,
                    pois_gamma = gamma_log_j)
    for_ru <- c(for_ru, list(phi_to_theta = phi_to_theta, log_j = log_j))
  }
  res <- do.call(rust::ru, for_ru)
  #
  # Sample from the conditional posterior distribution of the population
  # parameters given the hyperparameters and the data
  #
  temp <- switch(model,
                 binom_beta = binom_beta_cond_sim(res$sim_vals, data, n),
                 pois_gamma = pois_gamma_cond_sim(res$sim_vals, data, n))
  res <- c(res, temp)
  #
  # Add information about the model, the data and the prior
  #
  res$model <- model
  res$data <- data
  res$prior <- prior
  class(res) <- "hef"
  return(res)
}
