
<!-- README.md is generated from README.Rmd. Please edit that file -->
### Bayesian Analysis, No Gibbs

### What does bang do?

Please note that bang is under development.

Provides functions for the Bayesian analysis of some simple commonly-used models, without using Markov Chain Monte Carlo (MCMC) methods such as Gibbs sampling. The 'rust' package <https://cran.r-project.org/package=rust> is used to simulate a random sample from the required posterior distribution, using the ratio-of-uniforms method. Advantages of the ratio-of-uniforms method over MCMC in this context are that the user is not required to set tuning parameters nor to monitor convergence and a random posterior sample is produced. See the 'bang' website for more information, documentation and examples.

### A simple example

The `hef` function samples from the posterior distribution of the parameters of certain hierarchical exponential family models. The following code performs essentially the same analysis of the rat tumor data using a binomial-beta hierarchical model that appears in Section 5.3 of Gelman, A., Carlin, J. B., Stern, H. S. Dunson, D. B., Vehtari, A. and Rubin, D. B. (2013) Bayesian Data Analysis. Chapman & Hall / CRC. <https://www.stat.columbia.edu/~gelman/book>.

``` r
rat_res <- hef(model = "binom_beta", data = rat)
plot(rat_res)
```

### Installation

To install this development version from Github use:

``` r
install.packages("devtools")
library(devtools)
devtools::install_github("paulnorthrop/bang")
```
