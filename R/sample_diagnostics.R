# =========================== plot.hef ===========================

#' Plot diagnostics for a hef object
#'
#' \code{plot} method for class "hef".
#'
#' @param x an object of class "hef", a result of a call to \code{ru}.
#' @param y Not used.
#' @param type A character scalar.
#'
#'   If \code{type = "hyper"} then the posterior sample of \emph{all}
#'   hyperparameter values in \eqn{\phi} is plotted using
#'   \code{\link[rust]{plot.ru}}.
#'
#'   If \code{type = "ru"} then only the posterior sample generated using
#'   \code{\link[rust]{ru}} is plotted using \code{\link[rust]{plot.ru}}.
#'   This produces a different plot to \code{type = "hyper"} if \code{ru}
#'   is used only on a subset of \eqn{\phi}.  For example, this may be
#'   the case if \code{x} is the result of a call to \code{\link{hanova1}}.
#'   See vignette("revdbayes-anova-vignette", package = "bang") for
#'   information.
#'
#'   If \code{type = "pop_sim"} then sample(s) of the population-specific
#'   parameter vector \eqn{\theta} are plotted.  The population(s) included
#'   are determined by \code{which_pop}.  If \code{which_pop} has length
#'   one a histogram is produced using \code{\link[graphics]{hist}}.
#'   Otherwise, pairwise scatter plots are produced using
#'   \code{\link[graphics]{pairs}}
#'
#'   If \code{type = "pop_dens"} then posterior density estimates of
#'   population-specific parameter(s) are plotted. The population(s)
#'   included are determined by \code{which_pop}.
#'   The value of \code{super} determines whether or not these density
#'   estimates appear on one plot or in separate plots.
#' @param which_pop An integer vector.  If \code{type = "pop_sim"} or
#'   \code{type = "pop_dens"} then \code{which_pop} indicates which
#'   populations to include in the plot.
#' @param which_theta An integer vector indicating which component of
#'   of the population-specific parameter vector \eqn{\theta} to plot.
#'   Only relevant if \eqn{\theta} has length greater than one.
#' @param super A logical scalar.  If \code{type = "pop_dens"} then
#'   \code{super} indicates whether the density estimates should be
#'   superimposed in one graph or appear in separate graphs.
#' @param ... Additional arguments passed to \code{\link[rust]{plot.ru}},
#'   \code{\link[graphics]{hist}} or \code{\link[graphics]{pairs}}.
#'   In particular, \code{ru_scale = TRUE} produces a plot using the
#'   parameterisation used for ratio-of-uniforms sampling.
#' @examples
#' # Beta-binomial model, rat data
#' rat_res <- hef(model = "beta_binom", data = rat)
#'
#' # Hyperparameters alpha and beta
#' plot(rat_res)
#' # Parameterisation used for sampling
#' plot(rat_res, ru_scale = TRUE)
#'
#' # Population-specific posterior samples
#' plot(rat_res, type = "pop_sim")
#' plot(rat_res, type = "pop_sim", which_pop = c(1, 10))
#'
#' # Population-specific posterior densities
#' plot(rat_res, type = "pop_dens")
#' plot(rat_res, type = "pop_dens", which_pop = c(1, 10))
#' @seealso \code{\link[rust]{plot.ru}} for arguments that may be passed
#'   via ...., in particular \code{ru_scale}.
#' @export
plot.hef <- function(x, y, ..., type = c("hyper", "pop_sim", "pop_dens", "ru"),
                     which_pop = 1, which_theta = 1, super = TRUE) {
  if (!inherits(x, "hef")) {
    stop("use only with \"hef\" objects")
  }
  type <- match.arg(type)
  user_args <- list(...)
  # Use plot.ru() to plot the simulated hyperparameter values
  if (type == "hyper") {
    for_ru <- x
    if (is.null(user_args$ru_scale) || !user_args$ru_scale) {
      for_ru$d <- ncol(for_ru$sim_vals)
    }
    class(for_ru) <- "ru"
    plot(for_ru, ...)
  } else if (type == "ru" ) {
    for_ru <- x
    for_ru$sim_vals <- for_ru$sim_vals[, for_ru$ru]
    class(for_ru) <- "ru"
    plot(for_ru, ...)
  } else if (type == "pop_sim") {
    plot_data <- x$theta_sim_vals[, which_pop]
    var_names <- colnames(x$theta_sim_vals)[which_pop]
    if (length(which_pop) == 1) {
      if (is.null(user_args$xlab)) {
        graphics::hist(plot_data, prob = TRUE, main = "", xlab = "", ...)
        if (!is.null(var_names)) {
          graphics::title(xlab = parse(text = var_names))
        }
      } else {
        graphics::hist(plot_data, prob = TRUE, main = "", ...)
      }
    } else {
      if (is.null(user_args$labels)) {
        graphics::pairs(plot_data, labels = parse(text = var_names), ...)
      } else {
        graphics::pairs(plot_data, ...)
      }
    }
  } else {
    # Beta-binomial
    if (x$model == "beta_binom" ) {
      alpha <- x$sim_vals[, 1]
      beta <- x$sim_vals[, 2]
      # Numbers of successes
      yy <- x$data[, 1]
      # Number of trials
      nn <- x$data[, 2]
      #
      ep <- 1e-6
      temp <- seq(ep, 1 - ep, len = 100)
      pdfxi <- function(x, i) {
        mean(stats::dbeta(x, alpha + yy[i], beta + nn[i] - yy[i]))
      }
      density_values <- vapply(temp, pdfxi, 0, i = 1)
      plot(temp, density_values, ..., type = "l")
    }
#  } else {
#    # Gamma-Poisson
#    if (x$model == "gamma_pois") {
#      alpha <- x$sim_vals[, 1]
#      beta <- x$sim_vals[, 2]
#      y <- data[, 1]
#      off <- data[, 2]
#      i <- which_pop
#      temp <- range(x$theta_sim_vals[, i])
#      temp <- seq(temp[1], temp[2], len = 101)
#      stats::dgamma(temp, shape = alpha + y[i], rate = beta + off[i])
#    }
  }
}

# =========================== summary.hef ===========================

#' Summarizing hef objects
#'
#' \code{summary} method for class "hef".
#'
#' @param object an object of class "hef", a result of a call to \code{hef}.
#' @param ... Additional arguments passed on to \code{\link[rust]{summary.ru}}.
#' @param type A character scalar.
#'
#'   If \code{type = "hyper"} then the posterior samples of all hyperparameter
#'   values in \eqn{\phi} are summarised using, using
#'   \code{\link[rust]{summary.ru}}.
#'
#'   If \code{type = "pop"} then only posterior samples of the populations
#'   specifed in \code{which_pop} are summarised.
#'
#' @param which_pop An integer vector.  If \code{type = "pop"} then
#'   \code{which_pop} indicates which populations, i.e. which columns
#'   of \code{object$theta_sim_vals} to summarise, using
#'   \code{\link{summary}}.  The default is all populations.
#' @examples
#' # Beta-binomial model, rat data
#' rat_res <- hef(model = "beta_binom", data = rat)
#'
#' # Posterior summaries of the hyperparameters alpha and beta
#' summary(rat_res)
#'
#' # Posterior summaries of the binomial probability for rats 1 to 3
#' summary(rat_res, type = "pop", which_pop = 1:3)
summary.hef <- function(object, ..., type = c("hyper", "pop"),
                        which_pop = 1:ncol(object$theta_sim_vals)) {
  if (!inherits(object, "hef")) {
    stop("use only with \"hef\" objects")
  }
  type <- match.arg(type)
  if (type == "hyper") {
    for_ru <- object
    class(for_ru) <- "ru"
    posterior_summary <- summary(for_ru, ...)
    return(invisible(posterior_summary))
  } else {
    posterior_summary <- summary(object$theta_sim_vals[, which_pop,
                                                       drop = FALSE])
  }
  return(posterior_summary)
}

