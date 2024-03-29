% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/set_and_check_prior.R
\name{set_user_prior}
\alias{set_user_prior}
\title{Set a user-defined prior}
\usage{
set_user_prior(
  prior,
  ...,
  model = c("beta_binom", "gamma_pois", "anova1"),
  anova_d = 2
)
}
\arguments{
\item{prior}{An R function returning \strong{the log of} the prior density
for of (perhaps a subset of) the hyperparameter vector \eqn{\phi}.}

\item{...}{Further arguments giving the names and values of any
parameters involved in the function \code{prior}.}

\item{model}{A character string.  Abbreviated name of the model:
"beta_binom" for beta-binomial and "gamma_pois" for gamma-Poisson
(see \code{\link{hef}}), "anova1" for 1-way ANOVA (see
\code{\link{hanova1}}).}

\item{anova_d}{An integer scalar.  Only relevant if \code{model = anova1}.
If \code{anova_d = 2} then \code{prior} must return the log-prior
density for the standard deviations \eqn{(\sigma_\alpha, \sigma)}
and a normal prior with mean \code{mu0} and standard deviation
\code{sigma0} is used for \eqn{\mu}.  The values of \code{mu0 = 0} and
\code{sigma0 = Inf} are set in the call to \code{hanova1}, with default
values \code{mu0 = 0} and \code{sigma0 = Inf}.
If \code{anova_d = 3} then \code{prior} must return the log-prior
density for \eqn{(\mu, \sigma_\alpha, \sigma)}.}
}
\value{
A list of class \code{"bang_prior"}.  Will contain the component
  \code{prior}, the user-supplied function to evaluate the log of the prior,
  and any arguments supplied in ....
}
\description{
Constructs a user-defined prior distribution for use as the argument
\code{prior} in \code{\link{hef}} or \code{\link{hanova1}}.
}
\details{
For details of the hyperparameters in \eqn{\phi} see the
  \strong{Details} section of \code{\link{hef}} for the models
  \code{beta_binom} and \code{gamma_pois} and of \code{\link{hanova1}}
  for the model \code{anova1}.
}
\examples{
# User-defined prior, passing parameters
# (equivalent to prior = "gamma" with hpars = c(1, 0.01, 1, 0.01))
user_prior <- function(x, hpars) {
  return(dexp(x[1], hpars[1], log = TRUE) + dexp(x[2], hpars[2], log = TRUE))
}
user_prior_fn <- set_user_prior(user_prior, hpars = c(0.01, 0.01))
}
\seealso{
\code{\link{hef}} for hierarchical exponential family models.

\code{\link{hanova1}} for hierarchical one-way analysis of
  variance (ANOVA).
}
