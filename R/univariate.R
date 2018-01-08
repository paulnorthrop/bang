# Big question. Should the main argument to the prior function be a named list
# or a vector?  parm is a vector once it gets to the loglik function,
# but this is because I have manipulated the function this way.

# Ease for user, consistency with hef and hanova (and other packages?)
# Less likely to make mistakes when?

# Give opportunity to check the names of the parameters in the prior,
# using an argument to set_user_prior?
# "optional but advisible"  !!!!!!!!!!!!!
# par_names <- NULL
# Do check and then remove this before we get to logpost.

# Input lists
# 1. start
# 2. param (or "original" or "trans") - "original" is the default for densfun a function
# 3. as part of the prior

# Force user to set prior.  Take parameterisation from prior and check for
# consistency with densfun and start.

# 1. Return -Inf for out-of-bounds (and document this in iid())
# 2. Transformation: lambda = 0 for positive parameters
# 3. How does ru() use lower? My doc claims for all optimisations if rotate = FALSE
#    ... but my ru algorithm bit doesn't know about this!
#    ... and that's where we get the problems.
#  The solutions are 1. and/or 2 (belt and braces)

# logit transformation for probabilities?
# lambda vs trans = "user", i.e. set phi_to_theta and log_j like hef()

# Do param = "original" and param = "trans" for default distns
# ... but allow lambda to override, i.e. only use param if lambda is null
# NO.  Use a list of transformation

# dgamma: rate vs scale is tricky.  Define a new dgamma if gamma_rate = FALSE?

# dnbinom, size and mu or size and prob

# Set sensible default priors
# If prior is conjugate the simulate (more) directly and return pdf of posterior

# 1. no fixed_pars
# 2. lambda
# 3. fixed_pars

# Recognise dgamma as "gamma"? .. or spell out that initial estimates are not
# calculated if we send dgamma

# Hold parameter(s) fixed?
# ... has 3 functions: densfun, ru, fitdistr.  How to deal?
# Write my own init function in which parameters may be fixed
#   for *all* distributions (fitdistr only has this for the ones that use optim)

# ==================================== iid =================================== #
#
#' Bayesian Inference for Univariate Distributions
#'
#' Similar to fitdistr
#'
#' @param x A numeric vector of length at least one containing only
#'   \code{\link{finite}} values
#' @param densfun  Either a character string or a function returning a density
#'   evaluated at its first argument.
#'
#'   Distributions \code{"beta"}, \code{"cauchy"}, \code{"chi-squared"},
#'   \code{"exponential"}, \code{"f"}, \code{"gamma"}, \code{"geometric"},
#'   \code{"log-normal"}, \code{"lognormal"}, \code{"logistic"},
#'   \code{"negative binomial"}, \code{"normal"}, \code{"Poisson"},
#'   \code{"t"} and \code{"weibull"} are recognised, case being ignored.
#' @param start A named list giving the parameters to be optimized with
#'   initial values. This can be omitted for some of the named distributions
#'   and must be for others (see \strong{Details}).
#' @param nsim A numeric scalar. The number of values to be simulated from the
#'   posterior distribution.
#' @param prior Describe.
#' @param param Describe.
#' @param lambda A numeric vector.  rep_len
#' @param gamma_rate A logical scalar.
#' @param ... Additional parameters, either for \code{densfun} or
#'   \code{\link[rust]{ru}}.
#'   In particular, it can be used to specify bounds via lower
#'   or upper or both. If arguments of densfun (or the density function
#'   corresponding to a character-string specification) are included they will
#'   be held fixed.
#' @details Add details
#' @return Add return
#' @examples
#' x <- rgeom(10, 0.5)
#' # fitdistr(x, "geometric")
#' geom_prior <- set_user_prior(dbeta, shape1 = 1, shape2 = 1, log = TRUE,
#'                              model = "iid", par_names = "prob")
#' set.seed(47)
#' pjn1 <- iid(x, "geometric", prior = geom_prior)
#' set.seed(47)
#' pjn2 <- iid(x, "geometric", prior = "beta")
#' set.seed(47)
#' pjn3 <- iid(x, "geometric", prior = "default")
#' @seealso The \code{\link[rust]{ru}} function in the \code{\link{rust}}
#'   package for details of the arguments that can be passed to \code{ru}.
#' @seealso \code{\link{set_user_prior}} to set a user-defined prior.
#' @references Venables, W. N. & Ripley, B. D. (2002) Modern Applied
#'   Statistics with S.  Fourth Edition. Springer, New York.
#'   ISBN 0-387-95457-0. \url{http://www.stats.ox.ac.uk/pub/MASS4}.
#' @export
iid <- function(x, densfun, start, nsim = 1000, prior = "default",
                hpars = NULL, param = "trans",
                lambda = NULL, gamma_rate = TRUE, ...) {
  # Check that x and densfun are OK
  if (missing(x) || length(x) == 0L || mode(x) != "numeric") {
    stop("'x' must be a non-empty numeric vector")
  }
  if (any(!is.finite(x))) {
    stop("'x' contains missing or infinite values")
  }
  if (missing(densfun) || !(is.function(densfun) || is.character(densfun))) {
    stop("'densfun' must be supplied as a function or name")
  }
  if (missing(start)) {
    start <- NULL
  } else {
    if (!is.list(start)) {
      stop("'start' must be a named list")
    }
  }
  the_call <- match.call()
  n <- length(x)
  # Which arguments supplied in ... relate to densfun, optim and ru?
  ru_names <- names(formals(rust::ru))
  dots_list <- list(...)
  dots_names <- names(dots_list)
  not_ru_arg <- !is.element(dots_names, ru_names)
  fixed_pars <- dots_names[not_ru_arg]
  fixed_pars_list <- dots_list[not_ru_arg]
  ru_args <- dots_list[is.element(dots_names, ru_names)]
  # If densfun is a name then select the appropriate density function
  mydt <- function(x, m, s, df, log) {
    return(stats::dt((x - m) / s, df, log = TRUE) - log(s))
  }
  if (is.character(densfun)) {
    distname <- tolower(densfun)
    densfun <-
      switch(distname,
             "beta" = stats::dbeta,
             "cauchy" = stats::dcauchy,
             "chi-squared" = stats::dchisq,
             "exponential" = stats::dexp,
             "f" = stats::df,
             "gamma" = stats::dgamma,
             "geometric" = stats::dgeom,
             "log-normal" = stats::dlnorm,
             "lognormal" = stats::dlnorm,
             "logistic" = stats::dlogis,
             "negative binomial" = stats::dnbinom,
             "normal" = stats::dnorm,
             "poisson" = stats::dpois,
             "t" = mydt,
             "weibull" = stats::dweibull,
             NULL)
    if (is.null(densfun)) {
      stop("unsupported distribution")
    }
    if (is.null(start)) {
      start <- set_starting_values(distname, x, gamma_rate)
    }
  } else {
    distname <- "user-supplied"
  }
  # Remove from start any parameters to be kept fixed
  start <- start[!is.element(names(start), fixed_pars)]
  if (length(start) == 0L) {
    stop("All parameters are being held fixed")
  }
  nm <- names(start)
  print(prior)
  print(densfun)
  # Create a list that defines the prior and any parameters in the prior
  prior <- check_prior(prior, model = "iid", hpars = hpars,
                       distname = distname)
  print(prior)
  # Check that if the parameter names have been supplied in prior then
  # they match with those in start
  if (!is.null(prior$par_names)) {
    if (prior$par_names != nm) {
      stop("'prior$par_names' does not match 'names(start)'")
    }
    prior$par_names <- NULL
  }
#  print("nm")
#  print(nm)
#  print("start")
#  print(start)
#  print("prior")
#  print(prior$par_names)
  ## reorder arguments to densfun
  f <- formals(densfun)
  args <- names(f)
  m <- match(nm, args)
  if (any(is.na(m))) {
    stop("'start' includes names that are not arguments to 'densfun'")
  }
  # Set the function to calculate the log-likelihood
  formals(densfun) <- c(f[c(1, m)], f[-c(1, m)])
  dens <- function(parm, the_data, ...) densfun(the_data, parm, ...)
  if((l <- length(nm)) > 1L) {
    body(dens) <-
    parse(text = paste("densfun(the_data,",
                       paste("parm[", 1L:l, "]", collapse = ", "),
                       ", ...)"))
  }
  myfn <- function(parm, ...) sum(log(dens(parm, ...)))
  mylogfn <- function(parm, ...) sum(dens(parm, ..., log = TRUE))
  loglik_fn <- if ("log" %in% args) mylogfn else myfn
  d <- length(start)
  logpost <- function(parm, ...) {
    loglik <- suppressWarnings(loglik_fn(parm, ...))
    if (is.infinite(loglik) || is.nan(loglik)) {
      return(-Inf)
    }
    logprior <- do.call(prior$prior, c(list(parm), prior[-1]))
#    logprior <- 0
    return(loglik + logprior)
  }
  # Add to the list ru_args of arguments to pass to rust::ru()
  d <- length(start)
  if (is.null(ru_args$rotate) && d > 1) {
    ru_args$rotate <- TRUE
  }
  fr <- list(d = d, var_names = names(start))
  if (!is.null(lambda)) {
    ru_args$trans <- "BC"
    ru_args$lambda <- rep_len(lambda, d)
  }
  # Set ru_args$n_grid and ru_args$ep_bc to NULL just in case they have been
  # specified in ...
  ru_args$n_grid <- NULL
  ru_args$ep_bc <- NULL
  init <- as.numeric(start)
  for_ru <- c(list(logf = logpost, the_data = x), fixed_pars_list, fr,
              list(init = init, n = nsim), ru_args)
  res <- do.call(rust::ru, for_ru)
  res$call <- the_call
  class(res) <- "iid"
  return(res)
}

set_starting_values <- function(distname, x, gamma_rate) {
  if (distname %in% c("lognormal",  "log-normal")) {
    if (any(x <= 0)) {
      stop("Log Normal values must be positive")
    }
    n <- length(x)
    lx <- log(x)
    sd0 <- sqrt((n - 1) / n) * stats::sd(lx)
    mx <- mean(lx)
    return(list(meanlog = mx, sdlog = sd0))
  }
  if (distname == "normal") {
    n <- length(x)
    sd0 <- sqrt((n - 1) / n) * stats::sd(x)
    mx <- mean(x)
    return(list(mean = mx, sd = sd0))
  }
  if (distname == "poisson") {
    if (any(x < 0)) {
      stop("Poisson values must be non-negative integers")
    }
    return(list(lambda = mean(x)))
  }
  if (distname == "exponential") {
    if (any(x < 0)) {
      stop("Exponential values must be >= 0")
    }
    return(list(rate = 1 / mean(x)))
  }
  if (distname == "geometric") {
    if (any(x < 0)) {
      stop("Geometric values must be non-negative integers")
    }
    return(list(prob = 1 /(1 + mean(x))))
  }
  if (distname == "weibull" && is.null(start)) {
    ## log-Weibull is Gumbel, so start from that
    if (any(x <= 0)) {
      stop("Weibull values must be > 0")
    }
    lx <- log(x)
    m <- mean(lx)
    v <- stats::var(lx)
    shape <- 1.2 / sqrt(v)
    scale <- exp(m + 0.572 / shape)
    return(list(shape = shape, scale = scale))
  }
  if (distname == "gamma") {
    if (any(x < 0)) {
      stop("gamma values must be >= 0")
    }
    m <- mean(x)
    v <- stats::var(x)
    if (gamma_rate) {
      start <- list(shape = m ^ 2 / v, rate = m / v)
    } else {
      start <- list(shape = m ^ 2 / v, scale = v / m)
    }
    return(start)
  }
  if (distname == "negative binomial") {
    if (any(x < 0)) {
      stop("Negative binomial values must be non-negative integers")
    }
    m <- mean(x)
    v <- stats::var(x)
    size <- if (v > m) m ^ 2 /(v - m) else 100
    return(list(size = size, mu = m))
  }
  if (is.element(distname, c("cauchy", "logistic"))) {
    return(list(location = stats::median(x), scale = stats::IQR(x) / 2))
  }
  if (distname == "t") {
    return(list(m = stats::median(x), s = stats::IQR(x) / 2, df = 10))
  }
  return()
}

