#' bang: Bayesian Analysis with No Gibbs Sampling
#'
#' Performs Bayesian analyses using some simple commonly-used models.
#' The multivariate generalized ratio-of-uniforms method is used to simulate
#' random samples from the requied posterior distribution.
#' The user can either choose hyperparameter values of a default prior
#' distribution or specify their own prior distribution.
#'
#' @details Add some details here.
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
#' @name bang
#' @import methods
NULL

#' Rat tumor data
#'
#' Tumor incidence in 71 groups of rate from Tarone (1982).
#' The \code{rat} data frame has 71 rows and 2 columns.
#' Each row relates to a different group of rats.
#' The first column contains the number of rats with tumors.
#' The second column contains the total number of rats.
#'
#' @format A data frame with 71 rows and 2 columns.
#' @source Table 5.1 of Gelman, A., Carlin, J. B., Stern, H. S. Dunson, D. B.,
#'  Vehtari, A. and Rubin, D. B. (2013) \emph{Bayesian Data Analysis}.
#'  Chapman & Hall / CRC.
#'   \url{http://www.stat.columbia.edu/~gelman/book/data/rats.asc}
#' @references Tarone, R. E. (1982) The use of historical information in
#'   testing for a trend in proportions. \emph{Biometrics}, \strong{38},
#'   215-220. \url{https://doi.org/10.2307/2530304}
"rat"
