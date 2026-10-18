mk <- function() {
  m <- tsmod(marg("st", c(mu = 0, sigma = 0.01, df = 5)),
             sdvinecop(arfactcop(lambda = 0.9, w = 0.3), tail = pc("ivt", base = "joe"),
                       vt = c(0.6, 0.45), maxlag = 8, tautol = 1e-3))
  set.seed(11)
  list(m = m, x = tssim(m, 300))
}

test_that("fixed-parameter tsroll agrees with the model's own predictions and residuals", {
  s <- mk()
  fc <- tsroll(s$m, s$x, start = 250, levels = c(0.1, 0.5, 0.9), p = 0.05)
  expect_s4_class(fc, "tsforecast")
  expect_equal(fc@index, 250:300)
  # PIT equals the copula residual of the model
  r <- resid(s$m, s$x)
  expect_equal(fc@pit, r[250:300], tolerance = 1e-8)
  # quantiles equal tsmod predictions at the origin
  for (t in c(250L, 277L, 300L)) {
    q <- predict(s$m, s$x[seq_len(t - 1L)], type = "qf", x = c(0.1, 0.5, 0.9))
    expect_equal(unname(fc@quantile[t - 249L, ]), as.numeric(q), tolerance = 1e-6)
  }
  # ES is below VaR; VaR equals the 5% quantile prediction
  expect_true(all(fc@es < fc@var))
  rk <- tsrisk(s$m, s$x[1:249], p = 0.05)
  expect_equal(fc@var[1L, 1L], rk$var, tolerance = 1e-6)
  expect_equal(fc@es[1L, 1L], rk$es, tolerance = 1e-5)
})

test_that("the log predictive densities add up to the log-likelihood", {
  s <- mk()
  fc <- tsroll(s$m, s$x, start = 2, levels = 0.5, p = 0.05)
  u <- tm_pits(s$m@margin, s$x)
  # the copula density contributions plus margin log densities equal the joint log-likelihood
  cop_ll <- sum(fc@logdens - log(mdens(s$m@margin, s$x[-1L])))
  expect_equal(cop_ll, dv_loglik(u, dv_resolve(s$m@copula)), tolerance = 1e-8)
})

test_that("the generic per-origin pass matches the compiled pass", {
  s <- mk()
  a <- tsroll(s$m, s$x, start = 290, levels = c(0.2, 0.8), p = 0.05)
  op <- options(tscop.cpp = FALSE); on.exit(options(op))
  b <- tsroll(s$m, s$x, start = 290, levels = c(0.2, 0.8), p = 0.05)
  expect_equal(a@pit, b@pit, tolerance = 1e-6)
  expect_equal(a@quantile, b@quantile, tolerance = 1e-5)
})

test_that("refitting uses blocks, windows and keeps forecasts conditional on all data", {
  skip_on_cran()
  s <- mk()
  fc <- tsroll(s$m, s$x, start = 241, refit = 30, window = "rolling", width = 200, levels = 0.5, p = 0.05)
  expect_equal(nrow(fc@blocks), 2L)
  expect_equal(fc@blocks$window_from, c(41L, 71L))
  expect_equal(fc@blocks$window_to, c(240L, 270L))
  expect_equal(fc@blocks$status[1L], "given")
  expect_true(fc@blocks$status[2L] %in% c("ok", "failed"))
  expect_equal(length(fc@pit), 60L)
  expect_true(all(is.finite(fc@logdens)))
  expect_equal(nrow(fc@estimates), 2L)
})

test_that("tsroll validates its arguments", {
  s <- mk()
  expect_error(tsroll(s$m, s$x, start = 1), "start")
  expect_error(tsroll(s$m, s$x, start = 100, p = 1.5), "'p'")
  expect_error(tsroll(s$m, s$x, start = 100, refit = 0), "refit")
  expect_error(tsroll(s$m, s$x, start = 100, levels = c(0, 0.5)), "levels")
  expect_error(tsroll(s$m), "data")
})
