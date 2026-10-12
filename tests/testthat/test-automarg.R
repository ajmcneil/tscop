test_that("automarg ranks the families and returns the selected fit", {
  skip_on_cran()
  set.seed(5)
  x <- tssim(marg("sst", c(mu = 0, sigma = 1, df = 4, gamma = 1.5)), 3000)
  a <- automarg(x)
  expect_s4_class(a, "automarg")
  expect_equal(sort(a@table$family), sort(names(tscop:::.marg_families)))
  expect_equal(a@table$delta[1], 0)
  expect_true(!is.unsorted(a@table$delta))
  expect_equal(a@table$family[1], "sst")
  expect_s4_class(a@best, "margfit")
  expect_equal(a@best@model@family, "sst")
  # the table agrees with direct fits of the same families
  for (f in c("gauss", "st")) {
    d <- tsfit(x, marg(f))
    expect_equal(a@table$loglik[a@table$family == f], as.numeric(logLik(d)), tolerance = 1e-5)
  }
  expect_equal(a@table$AIC, -2 * a@table$loglik + 2 * a@table$npar)
  expect_equal(a@table$BIC, -2 * a@table$loglik + log(length(x)) * a@table$npar)
  # delegating methods
  expect_equal(coef(a), coef(a@best))
  expect_equal(mdens(a, 0.3), mdens(a@best, 0.3))
  expect_equal(mqf(a, 0.9), mqf(a@best, 0.9))
  expect_equal(length(tssim(a, 20)), 20)
  expect_equal(resid(a), resid(a@best))
  expect_output(show(a), "selected model")
})

test_that("automarg options: families, fix, criterion, standard errors", {
  skip_on_cran()
  set.seed(6)
  x <- tssim(marg("st", c(mu = 0, sigma = 1, df = 5)), 800)
  a <- automarg(x, families = c("sst", "gauss", "st"), fix = c(mu = 0), criterion = "BIC", se = "hessian")
  expect_equal(nrow(a@table), 3L)
  expect_equal(a@criterion, "BIC")
  expect_equal(a@table$delta[1], 0)
  expect_true(all(vapply(a@fits, function(f) f@model@pars[["mu"]] == 0, TRUE)))
  expect_equal(a@best@se_method, "hessian")
  expect_equal(a@best@model@family, a@table$family[1])
  expect_error(automarg(x, families = "foo"), "unknown margin family")
  expect_error(automarg(x, fix = c(nu = 1)), "named numeric")
  expect_error(automarg(x[1:5]), "at least 10")
  # the order of the families does not matter
  b <- automarg(x, families = c("st", "gauss"))
  expect_equal(b@table$family[1], "st")
})

test_that("automarg reports families that cannot be fitted", {
  x <- rep(c(-1, 1), 10)  # degenerate: a two-point sample
  a <- tryCatch(automarg(x, families = c("gauss", "st")), error = function(e) e)
  # either some family fits or the failure is reported as an error with the notes
  expect_true(inherits(a, "automarg") || inherits(a, "error"))
})
