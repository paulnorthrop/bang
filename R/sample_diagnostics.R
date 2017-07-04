# =========================== plot.hef ===========================

#' Plot diagnostics for a hef object
#'
#' \code{plot} method for class "hef".  Calls \code{\link[rust]{plot.ru}}.
#'
#' @param x an object of class "hef", a result of a call to \code{ru}.
#' @param y Not used.
#' @param ... Additional arguments passed on to \code{\link[rust]{plot.ru}}.
#' @examples
#' # Binomial-beta model, rat data
#' rat_res <- hef(model = "binom_beta", data = rat)
#' plot(rat_res)
#' plot(rat_res, ru_scale = TRUE)
#' @seealso \code{\link[rust]{plot.ru}} for the arguments that may be passed
#'   via ....
#' @export
plot.hef <- function(x, y, ...) {
  if (!inherits(x, "hef")) {
    stop("use only with \"hef\" objects")
  }
  for_ru <- x
  class(for_ru) <- "ru"
  plot(for_ru)
}

