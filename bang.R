#' bangs: Bayesian Analysis with No Gibbs Sampling
#'
#' Performs Bayesian analyses using some simple commonly-used models.
#' The multivariate generalized ratio-of-uniforms method is used to simulate
#' random samples from the requied posterior distribution.
#' The user can either choose hyperparameter values of a default prior
#' distribution or specify their own prior distribution.
#'
#' @references Northrop, P. J. (2016). rust: Ratio-of-Uniforms Simulation with
#'   Transformation. R package version 1.2.2.
#'   \url{https://cran.r-project.org/package=rust}.
#'
#' @seealso The \code{\link[rust]{ru}} and \code{\link[rust]{ru_rcpp}}
#'   functions in the \code{\link{rust}} package for details of the arguments
#'   that can be passed to \code{ru} via \code{rpost} and for the form of the
#'   object (of class "evprior") returned from \code{rpost}, which has the same
#'   structure as an object (of class "ru") returned by \code{ru} and
#'   \code{ru_rcpp}.
#' @docType package
#' @name bangs
#' @import methods
NULL

