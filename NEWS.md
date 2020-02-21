# bang 1.0.1

## Bug fixes and minor improvements

* pkgdown documentation at [https://paulnorthrop.github.io/bang/](https://paulnorthrop.github.io/bang/)

* Added a check that the argument `hpars` to `hef()` and to `hanova1()` is valid.

* Corrected typos in the vignette Hierarchical 1-way Analysis of Variance.  In the 9th line of the Appendix the expression given for ndot was incorrect.  It should simply be the sum of the sample sizes in the individual groups.  In the subsection "Marginal posterior density of phi" in the Appendix (ci-bi/ai)^2 has been corrected to (ci-bi^2/ai).  This is a typo: the ultimate expressions are correct.

* The summary method for class "hef" is now set up according to Section 8.1 of the R FAQ at (https://cran.r-project.org/doc/FAQ/R-FAQ.html).

* Corrected an error that meant that the S3 plot method plot.hef() was not resetting graphical parameters correctly on exit.

