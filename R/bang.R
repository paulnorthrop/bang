#' bang: Bayesian Analysis, No Gibbs
#'
#' Performs Bayesian analyses using some simple commonly-used models.
#' The multivariate generalized ratio-of-uniforms method is used to simulate
#' random samples from the required posterior distribution.
#' The user can either choose hyperparameter values of a default prior
#' distribution or specify their own prior distribution.
#'
#' @details Currently three conjugate hierarchical models are available:
#' beta-binomial, gamma-Poisson and 1-way Analysis of Variance (ANOVA).
#' The function \code{\link{hef}} produces random posterior samples from for the
#' beta-binomial and gamma-Poisson models.  The function \code{\link{hanova1}}
#' does this for the 1-way Analysis of Variance (ANOVA).
#' The \href{https://cran.r-project.org/package=rust}{rust} package is used to
#' produce these samples.
#'
#' See
#' vignette("bang-a-vignette", package = "bang") for a brief introduction
#' to the package and
#' vignette("bang-b-hef-vignette", package = "bang") and
#' vignette("bang-c-anova-vignette", package = "bang") for illustrations
#' of the use of the \code{\link{hef}} and \code{\link{hanova1}} functions.
#'
#' @seealso \code{\link{hef}} for hierarchical exponential family models.
#' @seealso \code{\link{hanova1}} for hierarchical one-way analysis of
#'   variance (ANOVA).
#' @seealso \code{\link{set_user_prior}} to set a user-defined prior.
#' @references Northrop, P. J. (2017). rust: Ratio-of-Uniforms Simulation with
#'   Transformation. R package version 1.2.3.
#'   \url{https://cran.r-project.org/package=rust}.
#'
#' @docType package
#' @name bang
#' @importFrom graphics plot
#' @importFrom bayesplot ppc_stat
#' @importFrom bayesplot ppc_stat_2d
#' @importFrom bayesplot ppc_dens_overlay
#' @importFrom bayesplot ppc_ecdf_overlay
#' @importFrom bayesplot ppc_dens
#' @importFrom bayesplot ppc_hist
#' @importFrom bayesplot ppc_boxplot
#' @importFrom bayesplot ppc_intervals
#' @importFrom bayesplot ppc_ribbon
#' @importFrom bayesplot pp_check
#' @importFrom bayesplot mcmc_areas
#' @importFrom bayesplot mcmc_intervals
#' @importFrom bayesplot mcmc_dens
#' @importFrom bayesplot mcmc_hist
NULL

#' Rat tumor data
#'
#' Tumor incidence in 71 groups of rate from Tarone (1982).
#' The matrix \code{rat} has 71 rows and 2 columns.
#' Each row relates to a different group of rats.
#' The first column (\code{y}) contains the number of rats with tumors.
#' The second column (\code{n}) contains the total number of rats.
#'
#' @format A matrix with 71 rows and 2 columns.
#' @source Table 5.1 of Gelman, A., Carlin, J. B., Stern, H. S. Dunson, D. B.,
#'  Vehtari, A. and Rubin, D. B. (2014) \emph{Bayesian Data Analysis},
#'  Chapman & Hall / CRC.
#'   \url{http://www.stat.columbia.edu/~gelman/book/data/rats.asc}
#' @references Tarone, R. E. (1982) The use of historical information in
#'   testing for a trend in proportions. \emph{Biometrics}, \strong{38},
#'   215-220. \url{https://doi.org/10.2307/2530304}
"rat"

#' Pump-failure data
#'
#' Data on pump failures from Gaver, D. P. and O'Muircheartaigh, I. G. (1987).
#' The matrix \code{pump} has 10 rows and 2 columns.
#' Each row relates to a different pump system.
#' The first column contains the number of pump failures.
#' The second column contains the length of operating time, in
#' thousands of hours.
#'
#' @format A matrix 10 rows and 2 columns.
#' @source Table 3 of Gaver, D. P. and O'Muircheartaigh, I. G. (1987).
#'   See also Gelfand, A. E. and Smith, A. F. M. (1990).
#' @references Gaver, D. P. and O'Muircheartaigh, I. G. (1987)
#'   Robust Empirical Bayes Analyses of Event Rates. \emph{Technometrics},
#'   \strong{29}, 1-15. \url{https://doi.org/10.1080/00401706.1987.10488178}
#' @references Gelfand, A. E. and Smith, A. F. M. (1990)
#' Sampling-Based Approaches to Calculating Marginal Densities.
#' \emph{Journal of the American Statistical Association}, \strong{85}(410),
#' 398-409. \url{https://doi.org/10.1080/01621459.1990.10476213}
#'
"pump"

#' Coagulation time data
#'
#' Coagulation time in seconds for blood drawn from 24 animals randomly
#' allocated to four different diets from Box, Hunter, and Hunter (1978).
#' The data frame \code{coagulation} has 24 rows and 2 columns.
#' Each row relates to a different animal.
#' Column 1 contains the coagulation times.
#' Column 2 contains a label for the type of diet: one of A, B, C or D.
#' @format A data frame with 24 rows and 2 columns.
#' @source Table 11.2 of Gelman, A., Carlin, J. B., Stern, H. S. Dunson, D. B.,
#'  Vehtari, A. and Rubin, D. B. (2014) \emph{Bayesian Data Analysis}.
#'  Chapman & Hall / CRC.
#'   \url{http://www.stat.columbia.edu/~gelman/book/}
#' @references Box, G. E. P., Hunter, W. G., and Hunter, J. S. (1978).
#'   Statistics for Experimenters. New York: Wiley.
"coagulation"

#' Mid 21st Century Global Temperature Projection Data
#'
#' Indices of global temperature change from late 20th century (1970-1999)
#' to mid 21st century (2020-2049) based on data produced by the Fifth
#' Coupled Model Intercomparison Project (CMIP5).
#'
#' The data frame \code{temp1} data frame has 270 rows and 4 columns.
#' Each row relates to a climate projection run from one of 38 different
#' General Circulation Models (GCMs) under a particular
#' Representative Concentration Pathway (RCP).
#' Use \code{table(temp1[, c("GCM", "RCP")])} to see the numbers of
#' runs under each RCP for each GCM.
#' See Van Vuuren et al (2011) for an overview of RCPs
#' and Northrop and Chandler (2014) for analyses of a similar
#' older dataset (CMIP3).
#' Column 1 contains the anomaly of the mean global temperature over
#' the time period 2020-2049 relative to the mean global temperature
#' over 1970-1999, i.e. the latter subtracted from the former.
#' Column 2 contains an abbreviation for the name of the climate modelling
#' research group and the GCM.
#' Column 3 contains the RCP in the format \code{rcpxx} where \code{xx}
#' is a radiative forcing level resulting from an anticipated future
#' greenhouse gas emissions.
#' Column 4 is the simulation run number.
#' @format A data frame with 270 rows and 4 columns.
#'   \itemize{
#'     \item{Column 1, \code{index}: }{anomaly of 2020-2049 mean relative to
#'       the 1970-1999 mean.}
#'     \item{Column 2, \code{GCM}: }{Abbreviated name of General Circulation
#'       Model.}
#'     \item{Column 3, \code{RCP}: }{Representative Concentration Pathway.
#'       One of rcp26, rcp45, rcp60, rcp85.}
#'     \item{Column 4, \code{run}: }{Simulation run number.}
#'  }
#' @source The raw data from which the indices are calculated are monthly
#'   CMIP5 scenario runs for global surface air temperature (tas)
#'   downloaded from the KNMI Climate Explorer (\url{https://climexp.knmi.nl/})
#'   on 4/3/2015.
#' @references Northrop, P.J. and R.E. Chandler (2014). Quantifying
#'   Sources of Uncertainty in Projections of Future Climate.
#'   \emph{Journal of Climate}, \strong{27}, 8793-8808.
#'   \url{https://doi.org/10.1175/JCLI-D-14-00265.1}
#' @references Van Vuuren, D. P., Edmonds, J., Kainuma, M., Riahi, K.
#'   Thomson, A., Hibbard, K., Hurtt, G. C., Kram, T., Krey, V.,
#'   Lamarque, J.-F. (2011). The representative concentration pathways:
#'   an overview. \emph{Climatic change}, \strong{109}, 5-31.
#'   \url{https://doi.org/10.1007/s10584-011-0148-z}
"temp1"

#' Late 21st Century Global Temperature Projection Data
#'
#' Indices of global temperature change from late 20th century (1970-1999)
#' to late 21st century (2069-2098) based on data produced by the Fifth
#' Coupled Model Intercomparison Project (CMIP5).
#'
#' The data frame \code{temp2} data frame has 270 rows and 4 columns.
#' Each row relates to a climate projection run from one of 38 different
#' General Circulation Models (GCMs) under a particular
#' Representative Concentration Pathway (RCP).
#' Use \code{table(temp2[, c("GCM", "RCP")])} to see the numbers of
#' runs under each RCP for each GCM.
#' See Van Vuuren et al (2011) for an overview of RCPs
#' and Northrop and Chandler (2014) for analyses of a similar
#' older dataset (CMIP3).
#' Column 1 contains the anomaly of the mean global temperature over
#' the time period 2069-2098 relative to the mean global temperature
#' over 1970-1999, i.e. the latter subtracted from the former.
#' Column 2 contains an abbreviation for the name of the climate modelling
#' research group and the GCM.
#' Column 3 contains the RCP in the format \code{rcpxx} where \code{xx}
#' is a radiative forcing level resulting from an anticipated future
#' greenhouse gas emissions.
#' Column 4 is the simulation run number.
#' @format A data frame with 270 rows and 4 columns.
#'   \itemize{
#'     \item{Column 1, \code{index}: }{anomaly of 2069-2098 mean relative to
#'       the 1970-1999 mean.}
#'     \item{Column 2, \code{GCM}: }{Abbreviated name of General Circulation
#'       Model.}
#'     \item{Column 3, \code{RCP}: }{Representative Concentration Pathway.
#'       One of rcp26, rcp45, rcp60, rcp85.}
#'     \item{Column 4, \code{run}: }{Simulation run number.}
#'  }
#' @source The raw data from which the indices are calculated are monthly
#'   CMIP5 scenario runs for global surface air temperature (tas)
#'   downloaded from the KNMI Climate Explorer (\url{https://climexp.knmi.nl/})
#'   on 4/3/2015.
#' @references Northrop, P.J. and R.E. Chandler (2014). Quantifying
#'   Sources of Uncertainty in Projections of Future Climate.
#'   \emph{Journal of Climate}, \strong{27}, 8793-8808.
#'   \url{https://doi.org/10.1175/JCLI-D-14-00265.1}
#' @references Van Vuuren, D. P., Edmonds, J., Kainuma, M., Riahi, K.
#'   Thomson, A., Hibbard, K., Hurtt, G. C., Kram, T., Krey, V.,
#'   Lamarque, J.-F. (2011). The representative concentration pathways:
#'   an overview. \emph{Climatic change}, \strong{109}, 5-31.
#'   \url{https://doi.org/10.1007/s10584-011-0148-z}
"temp2"

#' Weight Gained by Rats
#'
#' Data from an experiment to study weight gained by 10 rats fed on four different
#' diets, defined by a combination of the amount of protein (low and high)
#' and by the source of protein (beef and cereal).
#'
#' @format A data frame with 40 rows and 3 columns.
#'   \itemize{
#'     \item{Column 1, \code{source}: }{source of protein, a factor with levels
#'       \code{Beef} and \code{Cereal}.}
#'     \item{Column 2, \code{type}: }{amount of protein, a factor with levels
#'       \code{High} and \code{Low}.}
#'     \item{Column 3, \code{weightgain}: }{weight gained, in grams.}
#'  }
#' @source D. J. Hand, A. D. Lunn, K. J. McConway, and E. Ostrowski (1994).
#'   \emph{A Handbook of Small Datasets}, Chapman and Hall/CRC, London.
"weight_gain"
