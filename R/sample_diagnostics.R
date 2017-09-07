# =========================== plot.hef ===========================

#' Plot diagnostics for a hef object
#'
#' \code{plot} method for class "hef".  Calls \code{\link[rust]{plot.ru}}.
#'
#' @param x an object of class "hef", a result of a call to \code{ru}.
#' @param y Not used.
#' @param type A character scalar.
#'
#'   If \code{type = "hyper"} then the posterior sample of hyperparameter
#'   values is plotted using \code{\link[rust]{plot.ru}}.
#'
#'   If \code{type = "ru"} then only the posterior sample generated using
#'   \code{\link[rust]{ru}} is plotted using \code{\link[rust]{plot.ru}}.
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
#'   \code{type = "pop_dens"} \code{which_pop} indicates which populations
#'   to include in the plot.
#' @param which_theta An integer vector indicating which component of
#'   of the population-specific parameter vector \eqn{\theta} to plot.
#'   Only relevant if \eqn{\theta} has length greater than one.
#' @param super A logical scalar.  If \code{type = "pop_dens"} then
#'   \code{super} indicates whether the density estimates should be
#'   superimposed in one graph or appear in separate graphs.
#' @param ... Additional arguments passed to \code{\link[rust]{plot.ru}},
#'   \code{\link[graphics]{hist}} or \code{\link[graphics]{pairs}}.
#' @examples
#' # Beta-binomial model, rat data
#' rat_res <- hef(model = "beta_binom", data = rat)
#' # Hyperparameters alpha and beta
#' plot(rat_res)
#' plot(rat_res, ru_scale = TRUE)
#' # Population-specific parameter samples
#' plot(rat_res, type = "pop_sim")
#' plot(rat_res, type = "pop_sim", which_pop = c(1, 10))
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
#' \code{summary} method for class "hef"
#'
#' @param object an object of class "hef", a result of a call to \code{hef}.
#' @param ... Additional arguments passed on to \code{\link[rust]{summary.ru}}.
#' @examples
#' # Beta-binomial model, rat data
#' rat_res <- hef(model = "beta_binom", data = rat)
#' summary(rat_res)
summary.hef <- function(object, ...) {
  if (!inherits(object, "hef")) {
    stop("use only with \"hef\" objects")
  }
  for_ru <- object
  class(for_ru) <- "ru"
  summary(for_ru, ...)
}

