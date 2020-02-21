# =========================== plot.hef ===========================

#' Plot diagnostics for a hef object
#'
#' \code{plot} method for class "hef".
#'
#' @param x an object of class "hef", a result of a call to
#'   \code{\link[rust]{ru}}.
#' @param y Not used.
#' @param ... Additional arguments passed to \code{\link[rust]{plot.ru}},
#'   \code{\link[graphics]{hist}} or \code{\link[graphics]{pairs}}.
#'   In particular, \code{ru_scale = TRUE} produces a plot using the
#'   parameterization used for ratio-of-uniforms sampling.
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
#'      See vignette("bang-c-anova-vignette", package = "bang") for
#'      information.}
#'    \item{"pop": }{posterior samples and/or densities of the
#'      population-specific parameter \eqn{\theta} are plotted.  The
#'      population(s) included are determined by \code{which_pop} and the
#'      type of plot is determined by \code{plot_type}.
#'      If \code{plot_type} is not supplied then it is set to \code{"dens"}.}
#'  }
#' @param which_pop An integer vector or character scalar.
#'   If \code{params = "pop"} then \code{which_pop} indicates which
#'   populations to include in the plot.  If \code{which_pop} is supplied then
#'   \code{params} is set to "pop".  If \code{which_pop = "all"}
#'   then all populations are included.  If there are many populations then
#'   this may fail if \code{plot_type = "pairs"} and/or
#'   \code{one_plot = FALSE}.
#'
#' @param plot_type A character scalar that determines the type of plot
#'   produced when \code{params = "pop"}.  If \code{plot_type} is supplied
#'   then \code{params} is set automatically to \code{"pop"}.
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
#' @param one_plot,add_legend,legend_position,legend_text Only relevant if
#'   \code{plot_type = "dens"}.  If \code{one_plot = TRUE} then the estimated
#'   marginal posterior densities are plotted in the same graph and if
#'   \code{add_legend = TRUE} then a legend is added to this plot using
#'   \code{\link[graphics]{legend}} in the position indicated by
#'   the character scalar \code{legend_position}.
#'   A character vector \code{legend_text} may be used to override the
#'   default legend text.
#' @param num A numeric scalar.  If \code{plot_type == "dens"} or
#'   \code{plot_type == "both"} then \code{num} gives the number of
#'   points at which the marginal densities are evaluated to produce plots.
#' @section Examples:
#' See the examples in \code{\link{hef}} and \code{\link{hanova1}}.
#' @seealso \code{\link[rust]{plot.ru}} for arguments that may be passed
#'   via ...., in particular \code{ru_scale}.
#' @export
plot.hef <- function(x, y, ..., params = c("hyper", "ru", "pop"),
                     which_pop = NULL, plot_type = NULL, one_plot = FALSE,
                     add_legend = FALSE, legend_position = "topright",
                     legend_text = NULL, num = 100) {
  if (!inherits(x, "hef")) {
    stop("use only with \"hef\" objects")
  }
  params <- match.arg(params)
  if (params == "pop" & is.null(plot_type)) {
    plot_type <- "dens"
  }
  all_pop <- 1:ncol(x$theta_sim_vals)
  if (is.character(which_pop)) {
    which_pop <- match.arg(which_pop, "all")
    params <- "pop"
    which_pop <- all_pop
    if (is.null(plot_type)) {
      plot_type <- "dens"
    }
  } else if (is.numeric(which_pop)) {
    if (!all(which_pop %in% all_pop)) {
      stop("which_pop must be a subset of ", min(all_pop), ":", max(all_pop))
    }
    params <- "pop"
    if (is.null(plot_type)) {
      plot_type <- "dens"
    }
  } else {
    which_pop <- 1
  }
  if (!is.null(plot_type)) {
    plot_type <- match.arg(plot_type, c("sim", "dens", "both", "pairs"))
    if (plot_type == "pairs" & length(which_pop) == 1) {
      stop("If plot_type = ''pairs'' then length(which_pop) must be > 1")
    }
    if (one_plot & plot_type != "dens") {
      stop("one_plot = TRUE is not relevant unless plot_type = ''dens''")
    }
    params <- "pop"
  }
  # Save par settings so that we can reset them on exit
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
  var_names <- colnames(x$theta_sim_vals)[which_pop]
  # Set xlab, ylab and main
  if (!is.null(user_args$xlab)) {
    my_xlab <- rep_len(user_args$xlab, n_pop)
  } else {
    if (one_plot) {
      len_name <- nchar(var_names[1])
      my_xlab <- parse(text = substr(var_names, 1, len_name - 3))
    } else {
      my_xlab <- parse(text = var_names)
    }
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
  if (!is.null(legend_text)) {
    my_legend <- rep_len(legend_text, n_pop)
  } else {
    my_legend <- parse(text = var_names)
  }
  if (!is.null(user_args$lty)) {
    my_lty <- rep_len(user_args$lty, n_pop)
  } else {
    my_lty <- 1:n_pop
  }
  if (!is.null(user_args$col)) {
    my_col <- rep_len(user_args$col, n_pop)
  } else {
    my_col <- 1:n_pop
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
                        beta_binom = pde_beta_binom(x, which_pop, num),
                        gamma_pois = pde_gamma_pois(x, which_pop, num),
                        anova1 = pde_anova1(x, which_pop, num))
  }
  # Set the number of rows and columns in the plot
  rc <- n2mfrow(n_pop)
  oldpar <- graphics::par(mfrow = rc)
  on.exit(graphics::par(oldpar))
  pairwise_hist <- function(x, ..., xlab, ylab, main) {
    for (i in 1:n_pop) {
      graphics::hist(x[, i], prob = TRUE, main = my_main[i], xlab = my_xlab[i],
                     ylab = my_ylab[i], ...)
    }
  }
  pairwise_dens <- function(x, one_plot, add_legend, ..., xlab, ylab, main) {
    if (one_plot) {
      graphics::par(mfrow = c(1, 1))
      graphics::matplot(x$xx, x$yy, type = "l", main = my_main[1],
                        xlab = my_xlab[1], ylab = my_ylab[1], ...)
      fn_my_legend <- function(x, legend, ..., lty, col) {
        legend(x = x, legend = legend, lty = my_lty, col = my_col, ...)
      }
      if (add_legend) {
        fn_my_legend(x = legend_position, legend = my_legend, ...)
      }
    } else {
      for (i in 1:n_pop) {
        graphics::matplot(x$xx, x$yy[, i], type = "l", main = my_main[i],
                          xlab = my_xlab[i], ylab = my_ylab[i], ...)
      }
    }
  }
  pairwise_hist_dens <- function(x, y, ..., xlab, ylab, main) {
    for (i in 1:n_pop) {
      temp <- graphics::hist(x[, i], plot = FALSE)
      my_ylim <- c(0, max(temp$density, y$yy[, i]))
      graphics::hist(x[, i], prob = TRUE, main = my_main[i], xlab = my_xlab[i],
                     ylab = my_ylab[i], ylim = my_ylim, ...)
      graphics::lines(y$xx, y$yy[, i], ...)
    }
  }
  # Avoid leaving the function with daft value of pin in par if the plotting
  # fails because there are too many plots
  temp <- try(switch(plot_type,
         sim = pairwise_hist(plot_data, ...),
         dens = pairwise_dens(post_dens, one_plot = one_plot,
                              add_legend = add_legend, ...),
         both = pairwise_hist_dens(plot_data, post_dens, ...)),
         silent = FALSE)
  return(invisible())
}

# ================================= n2mfrow ===================================

n2mfrow <- function (nr.plots) {
  if (nr.plots <= 3)
    c(nr.plots, 1)
  else if (nr.plots <= 6)
    c((nr.plots + 1)%/%2, 2)
  else if (nr.plots <= 12)
    c((nr.plots + 2)%/%3, 3)
  else c(nrow <- ceiling(sqrt(nr.plots)), ceiling(nr.plots/nrow))
}

# =============================== pde_beta_binom ==============================

pde_beta_binom <- function(x, which_pop, num) {
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
  h_vals <- seq(min_p, max_p, len = num)
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

# =============================== pde_gamma_pois ==============================

pde_gamma_pois <- function(x, which_pop, num) {
  alpha <- x$sim_vals[, 1]
  beta <- x$sim_vals[, 2]
  # Counts
  y <- x$data[, 1]
  # Offset
  off <- x$data[, 2]
  # Set the (common) ranges of values on the horizontal axes
  ep <- 1e-6
  plot_data <- x$theta_sim_vals[, which_pop, drop = FALSE]
  min_val <- max(min(plot_data), ep)
  max_val <- max(plot_data)
  h_vals <- seq(min_val, max_val, len = num)
  # Function to estimate the posterior density for a given population
  pdfxi <- function(x, pop) {
    mean(stats::dgamma(x, shape = alpha + y[pop], rate = beta + off[pop]))
  }
  density_matrix <- matrix(NA, ncol = ncol(plot_data), nrow = length(h_vals))
  # Loop over which_pop
  for (i in 1:length(which_pop)) {
    density_matrix[, i] <- vapply(h_vals, pdfxi, 0, pop = which_pop[i])
  }
  return(list(xx = h_vals, yy = density_matrix))
}

# ================================= pde_anova1 ================================

pde_anova1 <- function(x, which_pop, num) {
  mu <- x$sim_vals[, 1]
  va <- x$sim_vals[, 2] ^ 2
  ve <- x$sim_vals[, 3] ^ 2
  # Responses
  resp <- x$data[, 1]
  # Explanatory factor
  fac <- x$data[, 2]
  # Set the (common) ranges of values on the horizontal axes
  plot_data <- x$theta_sim_vals[, which_pop, drop = FALSE]
  min_val <- min(plot_data)
  max_val <- max(plot_data)
  h_vals <- seq(min_val, max_val, len = num)
  # Function to estimate the posterior density for a given population
  ds <- x$summary_stats
  pdfxi <- function(x, pop) {
    ni <- ds$ni[pop]
    div <- va + ve / ni
    mean_alpha <- mu + va * (ds$ybari[pop] - mu) / div
    sd_alpha <- sqrt((va * ve / ni) / div)
    mean(stats::dnorm(x, mean = mean_alpha, sd = sd_alpha))
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
#' @param object an object of class "hef", a result of a call to
#'   \code{\link{hef}}.
#' @param ... Additional arguments passed on to \code{\link[rust]{summary.ru}}.
#' @param params A character scalar.
#'
#'   If \code{params = "hyper"} then the posterior samples of all
#'   hyperparameter values in \eqn{\phi} are summarized using
#'   \code{\link[rust]{summary.ru}}.
#'
#'   If \code{params = "pop"} then only posterior samples of the populations
#'   specified in \code{which_pop} are summarized.
#'
#' @param which_pop An integer vector.  If \code{params = "pop"} then
#'   \code{which_pop} indicates which populations, i.e. which columns
#'   of \code{object$theta_sim_vals} to summarize, using
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
#' @export
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
    class(posterior_summary) <- c("summary.hef", "summary.ru")
  } else {
    posterior_summary <- summary(object$theta_sim_vals[, which_pop,
                                                       drop = FALSE])
    class(posterior_summary) <- c("summary.hef", "table")
  }
  return(posterior_summary)
}

# ============================= print.summary.hef =============================

#' Print method for objects of class "summary.hef"
#'
#' \code{print} method for class "summary.hef".
#'
#' @param x an object of class "summary.hef", a result of a call to
#'   \code{\link{summary.hef}}.
#' @param ... Additional optional arguments to be passed to
#'   \code{\link{print}}.
#' @details What is printed depends on the argument \code{params} supplied
#'   to \code{\link{summary.hef}}.
#' @return The argument \code{x}, invisibly, as for all \code{\link{print}}
#'   methods.
#' @seealso \code{\link{summary.hef}}: \code{summary} method for class "hef".
#' @seealso \code{\link{hef}} for hierarchical exponential family models.
#' @seealso \code{\link{hanova1}} for hierarchical one-way analysis of
#'   variance (ANOVA).
#' @export
print.summary.hef <- function(x, ...) {
  if (!inherits(x, "summary.hef")) {
    stop("use only with \"summary.hef\" objects")
  }
  class(x) <- class(x)[2]
  print(x, ...)
  invisible(x)
}

# ================================= print.hef =================================

#' Print method for objects of class "hef"
#'
#' \code{print} method for class "hef".
#'
#' @param x an object of class "hef", a result of a call to
#'   \code{\link{hef}} or \code{\link{hanova1}}.
#' @param ... Additional optional arguments. At present no optional
#'   arguments are used.
#' @details Prints the original call to \code{\link{hef}} or
#'   \code{\link{hanova1}}, the name of the model and the number of populations
#'   in the hierarchical model.
#' @return The argument \code{x}, invisibly, as for all \code{\link{print}}
#'   methods.
#' @seealso \code{\link{hef}} for hierarchical exponential family models.
#' @seealso \code{\link{hanova1}} for hierarchical one-way analysis of
#'   variance (ANOVA).
#' @export
print.hef <- function(x, ...) {
  if (!inherits(x, "hef")) {
    stop("use only with \"hef\" objects")
  }
  cat("\n", "Call:", paste(deparse(x$call)), "\n")
  cat("Model:", x$model, "\n")
  cat("Number of populations:", ncol(x$theta_sim_vals), "\n")
  cat("\n")
  return(invisible(x))
}
