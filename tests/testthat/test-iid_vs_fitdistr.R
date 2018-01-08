context("iid, check that MAP equals MLE, if a flat prior is used")

my_seed <- 47
my_tol <- 1e-5
my_n <- 1

# Is MASS available?
got_MASS <- requireNamespace("MASS", quietly = TRUE)

set.seed(my_seed)

if(got_MASS) {
  # Geometric
  x <- rgeom(11, 0.5)
  geom_prior <- set_user_prior(dbeta, shape1 = 1, shape2 = 1, log = TRUE,
                               model = "iid", par_names = "prob")
  res1 <- iid(x, "geometric", prior = geom_prior, n = my_n)$f_mode
  res2 <- fitdistr(x, "geometric")$estimate
  test_that("iid geometric: MAP vs MASS fitdistr() MLE", {
    testthat::expect_equivalent(res1, res2, tolerance = my_tol)
  })
}
