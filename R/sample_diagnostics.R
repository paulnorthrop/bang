# =========================== plot.hef ===========================

#' Plot diagnostics for a hef object
#'
#' \code{plot} method for class "hef".
#'
#' @param x an object of class "hef", a result of a call to \code{ru}.
#' @param y Not used.
#' @param ... Additional arguments passed to \code{\link[rust]{plot.ru}},
#'   \code{\link[graphics]{hist}} or \code{\link[graphics]{pairs}}.
#'   In particular, \code{ru_scale = TRUE} produces a plot using the
#'   parameterisation used for ratio-of-uniforms sampling.
#' @param params A character scalar that determines to which parameters the
#'   plots relate.
#'  \itemize{
#'    \item{"hyper": }{the posterior sample of \emph{all} hyperparameter
#'      values in \eqn{\phi} is plotted using \code{\link[rust]{plot.ru}}.}
#'    \item{"ru": }{then only the posterior sample generated using
#'      \code{\link[rust]{ru}} is plotted using \code{\link[rust]{plot.ru}}.
#'      This produces a different plot to \code{params = "hyper"} if \code{ru}
#'      is used only on a subset of \eqn{\phi}.  For example, this may be
#'      the case if \code{x} is the result of a call to \code{\link{hanova1}}.
#'      See vignette("revdbayes-anova-vignette", package = "bang") for
#'      information.}
#'    \item{"pop": }{posterior samples and/or densities of the
#'      population-specific parameter \eqn{\theta} are plotted.  The
#'      population(s) included are determined by \code{which_pop} and the
#'      type of plot is determined by \code{plot_type}.}
#'  }
#' @param which_pop An integer vector.  If \code{params = "pop"} then
#'   \code{which_pop} indicates which populations to include in the plot.
#' @param plot_type A character scalar that determines the type of plot
#'   produced when \code{params = "pop"}.
#'  \itemize{
#'    \item{"sim": }{histograms of the posterior samples
#'      of \eqn{\theta} for the populations in \code{which_pop}.}
#'    \item{"dens": }{estimates of the marginal posterior
#'      densities of \eqn{\theta} for the populations in \code{which_pop}.}
#'    \item{"both": }{both the histograms and estimated posterior densities.}
#'    \item{"pairs": }{pairwise scatter plots of the posterior samples of
#'      \eqn{\theta} for the populations in \code{which_pop}, which must have
#'      length greater than one.}
#'  }
#' @param rows A numeric scalar.  If \code{plot_type} is either \code{"sim"},
#'   \code{"dens"} or \code{"both"} then this sets the number of rows of plots.
#'   If this is not supplied then it is set internally.
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
#' plot(rat_res, params = "pop", plot_type = "both")
#' plot(rat_res, params = "pop", plot_type = "both", which_pop = c(1, 2, 10))
#'
#' # Population-specific posterior densities
#' plot(rat_res, params = "pop")
#' plot(rat_res, params = "pop", which_pop = c(1, 10))
#' @seealso \code{\link[rust]{plot.ru}} for arguments that may be passed
#'   via ...., in particular \code{ru_scale}.
#' @export
plot.hef <- function(x, y, ..., params = c("hyper", "ru", "pop"),
                     which_pop = 1,
                     plot_type = c("sim", "dens", "both", "pairs"),
                     rows = NULL) {
  if (!inherits(x, "hef")) {
    stop("use only with \"hef\" objects")
  }
  params <- match.arg(params)
  plot_type <- match.arg(plot_type)
  if (plot_type == "pairs" & length(which_pop) == 1) {
    stop("If plot_type = ''pairs'' then length(which_pop) must be > 1")
  }
  user_args <- list(...)
  # Use plot.ru() to plot the simulated hyperparameter values
  if (params == "hyper") {
    for_ru <- x
    if (is.null(user_args$ru_scale) || !user_args$ru_scale) {
      for_ru$d <- ncol(for_ru$sim_vals)
    }
    class(for_ru) <- "ru"
    plot(for_ru, ...)
    return(invisible())
  }
  if (params == "ru" ) {
    for_ru <- x
    for_ru$sim_vals <- for_ru$sim_vals[, for_ru$ru]
    class(for_ru) <- "ru"
    plot(for_ru, ...)
    return(invisible())
  }
  # Otherwise, plot population values
  n_pop <- length(which_pop)
  plot_data <- x$theta_sim_vals[, which_pop, drop = FALSE]
  # Set xlab, ylab and main
  if (!is.null(user_args$xlab)) {
    my_xlab <- rep_len(user_args$xlab, n_pop)
  } else {
    var_names <- colnames(x$theta_sim_vals)[which_pop]
    my_xlab <- parse(text = var_names)
  }
  if (!is.null(user_args$ylab)) {
    my_ylab <- rep_len(user_args$ylab, n_pop)
  } else {
    my_ylab <- rep("density", n_pop)
  }
  if (!is.null(user_args$main)) {
    my_main <- rep_len(user_args$main, n_pop)
  } else {
    my_main = rep("", n_pop)
  }
  # Pairs
  if (plot_type == "pairs") {
    if (is.null(user_args$labels)) {
      graphics::pairs(plot_data, labels = my_xlab, ...)
    } else {
      graphics::pairs(plot_data, ...)
    }
    return(invisible())
  }
  # If we need estimates of marginal posterior densities
  if (plot_type == "dens" || plot_type == "both") {
    post_dens <- switch(x$model,
                        beta_binom = pde_beta_binom(x, which_pop))
  }
  # Set the number of rows and columns in the plot
  if (is.null(rows)) {
    rc <- n2mfrow(n_pop)
  }
  def.par <- graphics::par(no.readonly = TRUE)
  graphics::par(mfrow = rc)
  pairwise_hist <- function(x, ..., xlab, ylab, main) {
    for (i in 1:n_pop) {
      graphics::hist(x[, i], prob = TRUE, main = my_main[i], xlab = my_xlab[i],
                     ylab = my_ylab[i], ...)
    }
  }
  pairwise_dens <- function(x, ..., xlab, ylab, main) {
    for (i in 1:n_pop) {
      graphics::matplot(x$xx, x$yy[, i], type = "l", main = my_main[i],
                        xlab = my_xlab[i], ylab = my_ylab[i], ...)
    }
  }
  pairwise_hist_dens <- function(x, y, ..., xlab, ylab, main) {
    for (i in 1:n_pop) {
      temp <- graphics::hist(x[, i], plot = FALSE)
      my_ylim <- c(0, max(temp$density, y$yy[, i]))
      graphics::hist(x[, i], prob = TRUE, main = my_main[i], xlab = my_xlab[i],
                     ylab = my_ylab[i], ylim = my_ylim, ...)
      lines(y$xx, y$yy[, i], ...)
    }
  }
  # Avoid leaving the function with daft value of pin in par if the plotting
  # fails because there are too many plots
  temp <- try(switch(plot_type,
         sim = pairwise_hist(plot_data, ...),
         dens = pairwise_dens(post_dens, ...),
         both = pairwise_hist_dens(plot_data, post_dens, ...)),
         silent = FALSE)
  if (inherits(temp, "try-error")) {
    graphics::par(def.par)
  }
  return(invisible())
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

# ============================== n2mfrow ==============================

n2mfrow <- function (nr.plots) {
  if (nr.plots <= 3)
    c(nr.plots, 1)
  else if (nr.plots <= 6)
    c((nr.plots + 1)%/%2, 2)
  else if (nr.plots <= 12)
    c((nr.plots + 2)%/%3, 3)
  else c(nrow <- ceiling(sqrt(nr.plots)), ceiling(nr.plots/nrow))
}

# =========================== post_dens_est ===========================

pde_beta_binom <- function(x, which_pop) {
  alpha <- x$sim_vals[, 1]
  beta <- x$sim_vals[, 2]
  # Numbers of successes
  yy <- x$data[, 1]
  # Number of trials
  nn <- x$data[, 2]
  # Set the (common) ranges of values on the horizontal axes
  ep <- 1e-6
  plot_data <- x$theta_sim_vals[, which_pop, drop = FALSE]
  min_p <- max(min(plot_data), ep)
  max_p <- min(max(plot_data), 1 - ep)
  h_vals <- seq(min_p, max_p, len = 100)
  # Function to estimate the posterior density for a given population
  pdfxi <- function(x, pop) {
    mean(stats::dbeta(x, alpha + yy[pop], beta + nn[pop] - yy[pop]))
  }
  density_matrix <- matrix(NA, ncol = ncol(plot_data), nrow = length(h_vals))
  # Loop over which_pop
  for (i in 1:length(which_pop)) {
    density_matrix[, i] <- vapply(h_vals, pdfxi, 0, pop = which_pop[i])
  }
  return(list(xx = h_vals, yy = density_matrix))
}

# =========================== summary.hef ===========================

#' Summarizing hef objects
#'
#' \code{summary} method for class "hef".
#'
#' @param object an object of class "hef", a result of a call to \code{hef}.
#' @param ... Additional arguments passed on to \code{\link[rust]{summary.ru}}.
#' @param params A character scalar.
#'
#'   If \code{params = "hyper"} then the posterior samples of all hyperparameter
#'   values in \eqn{\phi} are summarised using, using
#'   \code{\link[rust]{summary.ru}}.
#'
#'   If \code{params = "pop"} then only posterior samples of the populations
#'   specifed in \code{which_pop} are summarised.
#'
#' @param which_pop An integer vector.  If \code{params = "pop"} then
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
#' summary(rat_res, params = "pop", which_pop = 1:3)
summary.hef <- function(object, ..., params = c("hyper", "pop"),
                        which_pop = 1:ncol(object$theta_sim_vals)) {
  if (!inherits(object, "hef")) {
    stop("use only with \"hef\" objects")
  }
  params <- match.arg(params)
  if (params == "hyper") {
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

