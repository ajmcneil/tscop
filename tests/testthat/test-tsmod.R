test_that("tsmod validates and shows", {
  m <- tsmod(marg("st", c(mu = 0, sigma = 1, df = 5)), arfactcop(lambda = 0.9, w = 0.6))
  expect_s4_class(m, "tsmod")
  expect_s4_class(tsmod("sst", arfactcop(order = 1))@margin, "marg")
  expect_error(tsmod(1, arfactcop(order = 1)), "margin")
  expect_error(tsmod(marg("st"), "foo"), "copula")
  expect_error(tsmod(marg("st"), marg("st")), "copula")
  expect_output(show(m), "tsmod")
  expect_equal(names(coef(m)), c("margin.mu", "margin.sigma", "margin.df", "copula.lambda1", "copula.w1"))
  expect_equal(unname(coef(tsmod(emarg(rnorm(20)), arfactcop(lambda = 0.5, w = 0.5)))), c(0.5, 0.5))
  # every copula class is accepted
  for (cop in list(armacop(ar = 0.5), arfactcop(order = 1), sdvinecop(arfactcop(order = 1)),
                   vtscop(arfactcop(order = 1), udp::vlinear(0.5)))) {
    expect_s4_class(tsmod(marg("gauss"), cop), "tsmod")
  }
})

test_that("tssim applies the margin quantile function to the copula series", {
  set.seed(1)
  m <- tsmod(marg("st", c(mu = 1, sigma = 2, df = 5)), arfactcop(lambda = 0.9, w = 0.6))
  x <- tssim(m, 3000)
  expect_length(x, 3000)
  expect_equal(median(x), 1, tolerance = 0.3)
  # the copula residuals are uniform and serially independent
  r <- resid(m, x)
  expect_gt(stats::ks.test(r, "punif")$p.value, 0.001)
  expect_lt(max(abs(acf(qnorm(r), lag.max = 10, plot = FALSE)$acf[-1])), 0.06)
  # margin PITs
  expect_equal(resid(m, x, type = "margin"), pt((x - 1) / 2, 5), tolerance = 1e-8)
  expect_equal(resid(m, x, scale = "normal"), qnorm(r))
  expect_equal(residuals(m, x), r)
  expect_error(tssim(tsmod(marg("st"), arfactcop(lambda = 0.9, w = 0.6)), 10), "unspecified")
  expect_equal(ktau(m, 5), ktau(arfactcop(lambda = 0.9, w = 0.6), 5))
})

test_that("two-stage estimation recovers margin and copula parameters", {
  skip_on_cran()
  set.seed(2)
  truth <- tsmod(marg("sst", c(mu = 0, sigma = 1, df = 6, gamma = 1.2)), arfactcop(lambda = 0.9, w = 0.6))
  x <- tssim(truth, 3000)
  f <- tsfit(x, tsmod("sst", arfactcop(order = 1)))
  expect_s4_class(f, "tsmodfit")
  expect_equal(names(coef(f)), c("margin.mu", "margin.sigma", "margin.df", "margin.gamma", "copula.lambda1", "copula.w1"))
  expect_equal(unname(coef(f)[c("copula.lambda1", "copula.w1")]), c(0.9, 0.6), tolerance = 0.1)
  expect_equal(unname(coef(f)[c("margin.sigma", "margin.gamma")]), c(1, 1.2), tolerance = 0.1)
  expect_equal(as.numeric(logLik(f)), as.numeric(logLik(f@margfit)) + as.numeric(logLik(f@copfit)))
  expect_equal(attr(logLik(f), "df"), 6L)
  expect_equal(f@nobs, 3000L)
  expect_equal(AIC(f), -2 * as.numeric(logLik(f)) + 2 * f@npar)
  expect_equal(f@model@margin@family, "sst")
  # methods of the fitted model
  expect_length(tssim(f, 50), 50)
  expect_equal(resid(f), resid(f@model, x))
  expect_gt(stats::ks.test(resid(f), "punif")$p.value, 0.001)
  expect_equal(ktau(f, 3), ktau(f@model, 3))
  expect_output(show(f), "tsmodfit")
  # a known margin: all parameters fixed
  fk <- tsfit(x, tsmod(marg("sst", c(mu = 0, sigma = 1, df = 6, gamma = 1.2), fix = c("mu", "sigma", "df", "gamma")),
                       arfactcop(order = 1)))
  expect_equal(names(coef(fk)), c("copula.lambda1", "copula.w1"))
  expect_equal(unname(coef(fk@model)), unname(coef(fk)))
})

test_that("semiparametric estimation uses the pseudo-observations of the empirical margin", {
  skip_on_cran()
  set.seed(3)
  truth <- tsmod(marg("st", c(mu = 0, sigma = 1, df = 4)), arfactcop(lambda = 0.85, w = 0.7))
  x <- tssim(truth, 2500)
  f <- tsfit(x, tsmod(emarg(x), arfactcop(order = 1)))
  expect_s4_class(f@margfit, "emarg")
  expect_equal(names(coef(f)), c("copula.lambda1", "copula.w1"))
  expect_equal(unname(coef(f)), c(0.85, 0.7), tolerance = 0.12)
  # the copula fit is that of the ranks
  g <- tsfit(rank(x) / (length(x) + 1), arfactcop(order = 1), pseudo = FALSE)
  expect_equal(unname(coef(f)), unname(coef(g)), tolerance = 1e-8)
  expect_equal(as.numeric(logLik(f)), as.numeric(logLik(g)), tolerance = 1e-8)
  expect_equal(f@npar, 2L)
  expect_output(show(f), "pseudo-likelihood")
  # simulation from the fitted semiparametric model uses empirical quantiles
  s <- tssim(f, 100)
  expect_true(all(s >= min(x) & s <= max(x)))
  # the bandwidth selector is kept when the margin is refitted
  expect_equal(tsfit(x, tsmod(emarg(rnorm(50), bw = "nrd0"), arfactcop(order = 1)))@margfit@bwsel, "nrd0")
})

test_that("with a Gaussian margin and a Gaussian copula the model is a Gaussian ARMA process", {
  skip_on_cran()
  set.seed(4)
  x <- 2 + 1.5 * as.numeric(arima.sim(list(ar = 0.6, ma = 0.3), 1500))
  f <- tsfit(x, tsmod("gauss", armacop(order = c(1, 1))))
  a <- arima(x, order = c(1, 0, 1))
  expect_equal(unname(coef(f)[c("copula.ar1", "copula.ma1")]), unname(coef(a)[c("ar1", "ma1")]), tolerance = 0.03)
  expect_equal(unname(coef(f)["margin.mu"]), unname(coef(a)["intercept"]), tolerance = 0.05)
  # arima maximises the same likelihood, so it can only be larger; two-stage is close to it
  expect_lte(as.numeric(logLik(f)), as.numeric(logLik(a)) + 1e-3)
  expect_gt(as.numeric(logLik(f)), as.numeric(logLik(a)) - 1)
})

test_that("tsmod works with the d-vine and v-transform copula processes and checks the data", {
  skip_on_cran()
  set.seed(5)
  mg <- marg("nig", c(mu = 0, sigma = 1, shape = 1.5, skew = 0.2))
  cop1 <- sdvinecop(arfactcop(lambda = 0.9, w = 0.6), pairs = list(pc("gumbel")))
  x1 <- tssim(tsmod(mg, cop1), 700)
  f1 <- tsfit(x1, tsmod(marg("nig"), sdvinecop(arfactcop(order = 1), pairs = list(pc("gumbel")))))
  expect_s4_class(f1@copfit, "sdvinecopfit")
  expect_equal(unname(coef(f1)[c("copula.lambda1", "copula.w1")]), c(0.9, 0.6), tolerance = 0.2)
  cop2 <- vtscop(arfactcop(lambda = 0.9, w = 0.6), udp::v2p(0.45, 0.8))
  x2 <- tssim(tsmod(marg("sst", c(mu = 0, sigma = 1, df = 8, gamma = 1)), cop2), 1500)
  f2 <- tsfit(x2, tsmod(marg("sst"), vtscop(arfactcop(order = 1), udp::v2p(0.5, 1))))
  expect_s4_class(f2@copfit, "vtscopfit")
  expect_true(all(c("copula.delta", "copula.kappa") %in% names(coef(f2))))
  expect_gt(stats::ks.test(resid(f2), "punif")$p.value, 0.001)
  expect_error(tsfit(c(x1[1:10], NA, x1[12:30]), tsmod(marg("nig"), arfactcop(order = 1))), "complete")
  expect_error(tsfit(x1[1:10], tsmod(marg("nig"), arfactcop(order = 1))), "at least 20")
  expect_error(tsfit("a", tsmod(marg("nig"), arfactcop(order = 1))), "numeric")
  expect_error(tsfit(x1, tsmod(marg("nig"), arfactcop(order = 1)), se = "bootstrap"), "not yet available")
})

test_that("ts attributes are restored in the residuals", {
  set.seed(6)
  x <- ts(tssim(tsmod(marg("gauss", c(mu = 0, sigma = 1)), arfactcop(lambda = 0.8, w = 0.5)), 300), start = 1990, frequency = 12)
  f <- tsfit(x, tsmod("gauss", arfactcop(order = 1)))
  expect_s3_class(resid(f), "ts")
  expect_equal(stats::tsp(resid(f)), stats::tsp(x))
})
