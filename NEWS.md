# bang 1.0.4.9000

## Bug fixes and minor improvements

* Corrected the description the `rat` data in its help file. Thank you to Léo Belzile for spotting this. [Issue 3](https://github.com/paulnorthrop/bang/issues/3#issue-2833715433)

# bang 1.0.4

## Bug fixes and minor improvements

* Fixed 2 \link{} targets with missing (rust) package anchors in Rd files for `hanova1()` and `hef()`.

# bang 1.0.3

## Bug fixes

* Fixed issues with the incorrect use of \itemize in some Rd files.

# bang 1.0.2

## Bug fixes and minor improvements

* In `beta_init_ests()` the initial estimate of `beta` has been corrected in the case `vp > mp * (1 - mp)`, where `vp` and `mp` are respectively the sample variance and sample mean of the input probabilities. Thank you to Thomas Richardson for spotting this.

* Create the help file for the package correctly, with alias bang-package.

* Activated 3rd edition of the `testthat` package

# bang 1.0.1

## Bug fixes and minor improvements

* pkgdown documentation at [https://paulnorthrop.github.io/bang/](https://paulnorthrop.github.io/bang/)

* Added a check that the argument `hpars` to `hef()` and to `hanova1()` is valid.

* Corrected typos in the vignette Hierarchical 1-way Analysis of Variance.  In the 9th line of the Appendix the expression given for ndot was incorrect.  It should simply be the sum of the sample sizes in the individual groups.  In the subsection "Marginal posterior density of phi" in the Appendix (ci-bi/ai)^2 has been corrected to (ci-bi^2/ai).  This is a typo: the ultimate expressions are correct.

* The summary method for class "hef" is now set up according to Section 8.1 of the R FAQ at (https://cran.r-project.org/doc/FAQ/R-FAQ.html).

* Corrected an error that meant that the S3 plot method plot.hef() was not resetting graphical parameters correctly on exit.

