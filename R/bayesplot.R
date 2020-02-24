# ========================= pp_check.hef ===========================

#' Posterior predictive checks for a hef object
#'
#' \code{pp_check} method for class "hef".  This provides an interface
#' to the functions that perform posterior predictive checks in the
#' \strong{bayesplot} package.  See \link[bayesplot]{PPC-overview} for
#' details of these functions.
#'
#' @aliases pp_check
#'
#' @param object An object of class "hef", a result of a call to
#'   \code{\link{hef}} or \code{\link{hanova1}}.
#' @param fun The plotting function to call.
#'   Can be any of the functions detailed at \link[bayesplot]{PPC-overview}.
#'   The "ppc_" prefix can optionally be dropped if fun is specified
#'   as a string.
#' @param raw Only relevant if \code{object$model = "beta_binom"} or
#'   \code{object$model = "gamma_pois"}.
#'   If \code{raw = TRUE} then the raw responses are used in the plots.
#'   Otherwise, the \emph{proportions} of successes are used in the
#'   \code{beta_binom} case and the \emph{exposure-adjusted rate} in the
#'   \code{gamma_pois} case.  In both cases the values used are
#'    \code{object$data[, 1] / object$data[, 2]} and the equivalent
#'    in \code{object$data_rep}.
#' @param nrep The number of predictive replicates to use.  If \code{nrep}
#'   is supplied then the first \code{nrep} rows of \code{object$data_rep}
#'   are used.  Otherwise, or if \code{nrep} is greater than
#'   \code{nrow(object$data_rep)}, then all rows are used.
#' @param ... Additional arguments passed on to bayesplot functions.
#'   See \strong{Examples} below.
#' @details For details of these functions see \link[bayesplot]{PPC-overview}.
#'   See also the vignettes
#'   \href{https://CRAN.R-project.org/package=bang}{Conjugate Hierarchical Models},
#'   \href{https://CRAN.R-project.org/package=bang}{Hierarchical 1-way Analysis of Variance}
#'   and the \strong{bayesplot} vignette
#'   \href{https://CRAN.R-project.org/package=bayesplot}{Graphical posterior predictive checks}.
#'
#'   The general idea is to compare the observed data \code{object$data}
#'   with a matrix \code{object$data_rep} in which each row is a
#'   replication of the observed data simulated from the posterior predictive
#'   distribution.  For greater detail see Chapter 6 of Gelman et al. (2013).
#'
#' @return A ggplot object that can be further customized using the
#'   \strong{ggplot2} package.
#' @seealso \code{\link{hef}} and \code{\link{hanova1}} for sampling
#'   from posterior distributions of hierarchical models.
#' @seealso \strong{bayesplot} functions \link[bayesplot]{PPC-overview},
#'   \link[bayesplot]{PPC-distributions},
#'   \link[bayesplot]{PPC-test-statistics},
#'   \link[bayesplot]{PPC-intervals},
#'   \link[bayesplot]{pp_check}.
#' @references Jonah Gabry (2016). bayesplot: Plotting for Bayesian
#' Models. R package version 1.1.0.
#' \url{https://CRAN.R-project.org/package=bayesplot}
#' @references Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B.,
#' Vehtari, A., and Rubin, D. B. (2013). \emph{Bayesian Data Analysis}.
#' Chapman & Hall/CRC Press, London, third edition. (Chapter 6).
#' \url{http://www.stat.columbia.edu/~gelman/book/}
#' @examples
#' ############################ Beta-binomial #################################
#'
#' # ------------------------- Rat tumor data ------------------------------- #
#'
#' rat_res <- hef(model = "beta_binom", data = rat, nrep = 50)
#'
#' # Overlaid density estimates
#' pp_check(rat_res)
#' \donttest{
#' # Overlaid distribution function estimates
#' pp_check(rat_res, fun = "ecdf_overlay")
#' }
#' # Multiple histograms
#' pp_check(rat_res, fun = "hist", nrep = 8)
#' \donttest{
#' # Multiple boxplots
#' pp_check(rat_res, fun = "boxplot")
#' # Predictive medians vs observed median
#' pp_check(rat_res, fun = "stat", stat = "median")
#' }
#' # Predictive (mean, sd) vs observed (mean, sd)
#' pp_check(rat_res, fun = "stat_2d", stat = c("mean", "sd"))
#'
#' ############################ Gamma-Poisson #################################
#'
#' # ------------------------ Pump failure data ------------------------------ #
#'
#' pump_res <- hef(model = "gamma_pois", data = pump, nrep = 50)
#'
#' \donttest{
#' # Overlaid density estimates
#' pp_check(pump_res)
#' # Predictive (mean, sd) vs observed (mean, sd)
#' pp_check(pump_res, fun = "stat_2d", stat = c("mean", "sd"))
#' }
#'
#' ###################### One-way Hierarchical ANOVA ##########################
#'
#' #----------------- Late 21st Century Global Temperature Data ------------- #
#'
#' RCP26_2 <- temp2[temp2$RCP == "rcp26", ]
#' temp_res <- hanova1(resp = RCP26_2[, 1], fac = RCP26_2[, 2], nrep = 50)
#' \donttest{
#' # Overlaid density estimates
#' pp_check(temp_res)
#' # Predictive (mean, sd) vs observed (mean, sd)
#' pp_check(temp_res, fun = "stat_2d", stat = c("mean", "sd"))
#' }
#' @export pp_check
#' @export
pp_check.hef <- function(object, fun = NULL, raw = FALSE, nrep = NULL, ...) {
  if (!inherits(object, "hef")) {
    stop("use only with \"hef\" objects")
  }
  if (is.null(object[["data_rep"]])) {
    stop("data_rep is NULL: call hef() or hanova1() again supplying nrep")
  }
  if (is.null(fun)) {
    fun <- "dens_overlay"
  }
  # Extract the observed data and the simulated replicates
  y <- object[["data"]]
  yrep <- object[["data_rep"]]
  # Keep only nrep replicates (or all replicates if nrep is too large)
  yrep <- yrep[1:min(nrep, nrow(yrep)), , drop = FALSE]
  # If the responses are binomial or Poisson then use the proportions of
  # successes (binomial) or the exposure-adjusted rate (Poisson) unless
  # the user tells us not to by setting raw = TRUE
  if (object$model %in% c("beta_binom", "gamma_pois") & !raw) {
    divisor <- y[, 2]
    y <- y[, 1] / divisor
    yrep <- t(t(yrep) / divisor)
  } else {
    y <- y[, 1]
  }
  # Call the bayesplot function ppc_fun
  return(bayesplot::pp_check(y, yrep, fun = fun, ...))
}
