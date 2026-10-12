models <- function() {
  list(
    gauss_arma = tsmod(marg("gauss", c(mu = 2, sigma = 1.5)), armacop(ar = 0.7, ma = 0.2)),
    st_afc = tsmod(marg("st", c(mu = 0, sigma = 1, df = 5)), arfactcop(lambda = 0.9, w = 0.6)),
    sst_dvine = tsmod(marg("sst", c(mu = 0, sigma = 1, df = 6, gamma = 1.2)),
                      sdvinecop(arfactcop(lambda = 0.9, w = 0.6), pairs = list(pc("gumbel")))),
    nig_vt = tsmod(marg("nig", c(mu = 0, sigma = 1, shape = 1.5, skew = 0.2)),
                   vtscop(arfactcop(lambda = 0.9, w = 0.6), udp::v2p(0.45, 0.8)))
  )
}

test_that("predictive df, qf and density are consistent for every copula class", {
  set.seed(1)
  for (nm in names(models())) {
    m <- models()[[nm]]
    x <- tssim(m, 120)
    ref <- stats::median(x)
    g <- ref + c(-3, -1, 0, 0.7, 2, 4) * sd(x)
    df <- predict(m, x, type = "df", x = g)
    expect_true(all(diff(df) > 0), info = nm)
    expect_equal(predict(m, x, type = "qf", x = df), g, tolerance = 1e-5, info = nm)
    e <- 1e-5
    num <- (predict(m, x, type = "df", x = g + e) - predict(m, x, type = "df", x = g - e)) / (2 * e)
    expect_equal(predict(m, x, type = "dens", x = g), num, tolerance = 1e-4, info = nm)
    expect_equal(integrate(function(v) predict(m, x, type = "dens", x = v), -Inf, Inf, rel.tol = 1e-6)$value, 1,
                 tolerance = 1e-4, info = nm)
  }
})

test_that("predictive PITs of simulated one-step predictions are uniform", {
  skip_on_cran()
  set.seed(2)
  m <- models()$st_afc
  X <- vapply(1:600, function(i) tssim(m, 40), numeric(40))
  pit <- vapply(1:600, function(i) predict(m, X[1:39, i], type = "df", x = X[40, i]), numeric(1))
  expect_gt(stats::ks.test(pit, "punif")$p.value, 0.001)
})

test_that("the full-model predictive functions are the margin and copula ones composed", {
  set.seed(3)
  m <- models()$sst_dvine
  x <- tssim(m, 100)
  u <- mdf(m@margin, x)
  g <- c(-1, 0.3, 2)
  expect_equal(predict(m, x, type = "df", x = g), predict(m@copula, u, type = "df", x = mdf(m@margin, g)))
  p <- c(0.05, 0.5, 0.95)
  expect_equal(predict(m, x, type = "qf", x = p), mqf(m@margin, predict(m@copula, u, type = "qf", x = p)))
  # a Gaussian ARMA model: the predictive distribution is normal with the known moments
  mg <- models()$gauss_arma
  xg <- tssim(mg, 80)
  mo <- predict(mg@copula, pnorm((xg - 2) / 1.5), type = "moments")
  expect_equal(predict(mg, xg, type = "df", x = 3), pnorm((3 - 2) / 1.5, mo$mean, mo$sd), tolerance = 1e-8)
  expect_equal(predict(mg, xg, h = 3, type = "qf", x = 0.3),
               2 + 1.5 * qnorm(0.3, predict(mg@copula, pnorm((xg - 2) / 1.5), h = 3, type = "moments")$mean,
                               predict(mg@copula, pnorm((xg - 2) / 1.5), h = 3, type = "moments")$sd), tolerance = 1e-8)
  expect_error(predict(m, x, type = "qf", x = 1.2), "probabilities")
  expect_error(predict(m, x, type = "df"), "'x'")
  expect_error(predict(m, c(x, NA), type = "df", x = 0), "complete")
  expect_error(predict(m, x, h = 2, type = "df", x = 0), "h = 1")
})

test_that("predict works for fitted models and empirical margins", {
  skip_on_cran()
  set.seed(4)
  truth <- tsmod(marg("st", c(mu = 0, sigma = 1, df = 5)), arfactcop(lambda = 0.9, w = 0.6))
  x <- tssim(truth, 1500)
  f <- tsfit(x, tsmod("st", arfactcop(order = 1)))
  expect_equal(predict(f, type = "qf", x = c(0.1, 0.9)), predict(f@model, x, type = "qf", x = c(0.1, 0.9)))
  fs <- tsfit(x, tsmod(emarg(x), arfactcop(order = 1)))
  q <- predict(fs, type = "qf", x = c(0.1, 0.5, 0.9))
  expect_true(all(diff(q) > 0))
  # the parametric and semiparametric predictive quantiles are close
  expect_equal(q, predict(f, type = "qf", x = c(0.1, 0.5, 0.9)), tolerance = 0.2)
  expect_equal(predict(fs, type = "df", x = q), c(0.1, 0.5, 0.9), tolerance = 0.02)
})

test_that("tsrisk gives the predictive value-at-risk and expected shortfall", {
  set.seed(5)
  m <- models()$sst_dvine
  x <- tssim(m, 150)
  r <- tsrisk(m, x, p = c(0.01, 0.05))
  expect_s3_class(r, "data.frame")
  expect_equal(names(r), c("p", "tail", "h", "var", "es"))
  expect_equal(r$var, predict(m, x, type = "qf", x = c(0.01, 0.05)))
  expect_true(all(r$es < r$var))
  # expected shortfall by direct integration of the predictive quantile function
  qf <- function(s) predict(m, x, type = "qf", x = s)
  es <- c(integrate(qf, 0, 0.01)$value / 0.01, integrate(qf, 0, 0.05)$value / 0.05)
  expect_equal(r$es, es, tolerance = 1e-4)
  # upper tail
  ru <- tsrisk(m, x, p = 0.05, tail = "upper")
  expect_equal(ru$var, predict(m, x, type = "qf", x = 0.95))
  expect_equal(ru$es, integrate(qf, 0.95, 1)$value / 0.05, tolerance = 1e-4)
  expect_gt(ru$es, ru$var)
  # Gaussian margin and copula: closed form
  mg <- models()$gauss_arma
  xg <- tssim(mg, 80)
  mo <- predict(mg@copula, pnorm((xg - 2) / 1.5), type = "moments")
  rg <- tsrisk(mg, xg, p = 0.05)
  expect_equal(rg$es, 2 + 1.5 * (mo$mean - mo$sd * dnorm(qnorm(0.05)) / 0.05), tolerance = 1e-6)
  expect_equal(rg$var, 2 + 1.5 * (mo$mean + mo$sd * qnorm(0.05)), tolerance = 1e-8)
  # fitted models, horizons and input checks
  f <- tsfit(xg, tsmod("gauss", armacop(order = c(1, 0))))
  expect_equal(tsrisk(f, p = 0.05)$var, tsrisk(f@model, xg, p = 0.05)$var)
  expect_equal(tsrisk(mg, xg, p = 0.05, h = 3)$h, 3)
  expect_error(tsrisk(mg, xg, p = 1.5), "probabilities")
  expect_error(tsrisk(marg("st"), xg), "tsmod")
})
