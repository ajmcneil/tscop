# ---- ktau -------------------------------------------------------------

test_that("ktau equals (2/pi) asin of the PACF from stats::ARMAacf", {
  cases <- list(
    list(ar = 0.6, ma = numeric(0)),
    list(ar = c(0.5, 0.3), ma = 0.4),
    list(ar = 0.95, ma = -0.85),
    list(ar = numeric(0), ma = c(0.4, -0.3)),
    arma_expand(ar = 0.4, ma = 0.3, sar = 0.5, sma = -0.4, period = 4)
  )
  for (cs in cases) {
    mod <- armacop(ar = cs$ar, ma = cs$ma)
    ref <- (2 / pi) * asin(as.numeric(ARMAacf(cs$ar, cs$ma, lag.max = 15, pacf = TRUE)))
    expect_equal(ktau(mod, 15), ref, tolerance = 1e-8)
  }
  m <- armacop(ar = 0.3, sma = 0.5, period = 4)
  e <- arma_coef_expanded(m)
  expect_equal(ktau(m, 10),
               (2 / pi) * asin(as.numeric(ARMAacf(e$ar, e$ma, lag.max = 10, pacf = TRUE))),
               tolerance = 1e-8)
})

test_that("ktau of an AR(p) vanishes beyond lag p and has length lagmax", {
  m <- armacop(ar = c(0.5, 0.2))
  k <- ktau(m, 6)
  expect_length(k, 6)
  expect_equal(k[3:6], rep(0, 4), tolerance = 1e-8)
  expect_equal(ktau(armacop(ar = 0.6), 1), (2 / pi) * asin(0.6))
  expect_length(ktau(m), 20)
})

test_that("ktau works on fitted objects and validates lagmax", {
  set.seed(1)
  f <- tsfit(tssim(armacop(ar = 0.5), 300), armacop(order = c(1, 0)))
  expect_equal(ktau(f, 4), ktau(f@model, 4))
  expect_error(ktau(armacop(ar = 0.5), 0), "lagmax")
  expect_error(ktau(armacop(order = c(1, 0))), "unspecified")
})

test_that("ktau agrees with tscopula::kendall", {
  skip_if_no_oracle()
  old <- tscopula::armacopula(list(ar = 0.5, ma = 0.4))
  expect_equal(ktau(armacop(ar = 0.5, ma = 0.4), 8), as.numeric(tscopula::kendall(old, 8)),
               tolerance = 1e-8)
})

# ---- resid ------------------------------------------------------------

test_that("residuals are the Cholesky-standardised normal scores", {
  set.seed(2)
  for (cs in list(list(ar = 0.5, ma = 0.4), list(ar = numeric(0), ma = 0.6),
                  list(ar = c(0.5, 0.3), ma = numeric(0)))) {
    mod <- armacop(ar = cs$ar, ma = cs$ma)
    u <- tssim(mod, 50)
    L <- t(chol(toeplitz(acf_ref(cs$ar, cs$ma, 49))))
    ref <- as.numeric(forwardsolve(L, qnorm(u)))
    expect_equal(resid(mod, u, scale = "normal"), ref, tolerance = 1e-8)
    expect_equal(resid(mod, u), pnorm(ref), tolerance = 1e-8)
    expect_equal(residuals(mod, u), resid(mod, u))
  }
})

test_that("residuals of data from the model are iid standard normal", {
  set.seed(3)
  mod <- armacop(ar = 0.7, ma = -0.3)
  u <- tssim(mod, 20000)
  r <- resid(mod, u, scale = "normal")
  expect_lt(abs(mean(r)), 0.03)
  expect_equal(var(r), 1, tolerance = 0.05)
  expect_lt(max(abs(acf(r, lag.max = 10, plot = FALSE)$acf[-1])), 0.03)
  ru <- resid(mod, u)
  expect_equal(mean(ru), 0.5, tolerance = 0.02)
})

test_that("residuals handle NA, fitted defaults and time series attributes", {
  set.seed(4)
  mod <- armacop(ar = 0.5)
  u <- tssim(mod, 60)
  L <- t(chol(toeplitz(acf_ref(0.5, numeric(0), 59))))
  ref <- as.numeric(forwardsolve(L, qnorm(u)))
  u2 <- u
  u2[30] <- NA
  r <- resid(mod, u2, scale = "normal")
  expect_true(is.na(r[30]))
  expect_equal(sum(is.na(r)), 1)
  expect_equal(r[1:29], ref[1:29], tolerance = 1e-8)
  # fitted object defaults to its own data
  f <- tsfit(u, armacop(order = c(1, 0)))
  expect_equal(resid(f), resid(f@model, u))
  # ts attributes
  ut <- ts(u, start = c(2000, 1), frequency = 12)
  rt <- resid(mod, ut)
  expect_s3_class(rt, "ts")
  expect_equal(tsp(rt), tsp(ut))
  expect_s3_class(resid(tsfit(ut, armacop(order = c(1, 0)))), "ts")
  skip_if_not_installed("zoo")
  uz <- zoo::zoo(u, order.by = as.Date("2020-01-01") + 0:59)
  rz <- resid(mod, uz)
  expect_s3_class(rz, "zoo")
  expect_equal(zoo::index(rz), zoo::index(uz))
})

test_that("residuals agree with tscopula once the filter has settled", {
  skip_if_no_oracle()
  set.seed(5)
  mod <- armacop(ar = 0.5, ma = 0.4)
  u <- tssim(mod, 80)
  old <- tscopula::armacopula(list(ar = 0.5, ma = 0.4))
  ro <- tscopula:::resid_armacopula(old, u)
  rn <- resid(mod, u, scale = "normal")
  expect_equal(ro[21:80], rn[21:80], tolerance = 1e-6)
})

# ---- predict ----------------------------------------------------------

test_that("moments match exact Gaussian conditioning, for several horizons", {
  set.seed(6)
  for (cs in list(list(ar = 0.5, ma = 0.4), list(ar = numeric(0), ma = c(0.5, 0.2)))) {
    mod <- armacop(ar = cs$ar, ma = cs$ma)
    u <- tssim(mod, 30)
    acf <- acf_ref(cs$ar, cs$ma, 40)
    pr <- predict(mod, u, h = c(1, 2, 5))
    expect_equal(pr$h, c(1, 2, 5))
    for (i in 1:3) {
      ref <- cond_moments(qnorm(u), acf, pr$h[i])
      expect_equal(c(pr$mean[i], pr$sd[i]), unname(ref), tolerance = 1e-8)
    }
  }
  # seasonal
  mod <- armacop(ar = 0.3, sma = 0.5, period = 4)
  e <- arma_coef_expanded(mod)
  u <- tssim(mod, 40)
  pr <- predict(mod, u, h = 3)
  ref <- cond_moments(qnorm(u), acf_ref(e$ar, e$ma, 60), 3)
  expect_equal(c(pr$mean, pr$sd), unname(ref), tolerance = 1e-8)
})

test_that("AR(1) predictions have the closed form", {
  mod <- armacop(ar = 0.6)
  u <- tssim(mod, 25)
  pr <- predict(mod, u, h = 1:4)
  xn <- qnorm(u[25])
  expect_equal(pr$mean, 0.6^(1:4) * xn, tolerance = 1e-8)
  expect_equal(pr$sd, sqrt(1 - 0.6^(2 * (1:4))), tolerance = 1e-8)
})

test_that("predictions condition correctly on trailing and interior NA", {
  set.seed(7)
  mod <- armacop(ar = 0.5, ma = 0.4)
  u <- tssim(mod, 30)
  u[c(10, 29, 30)] <- NA
  acf <- acf_ref(0.5, 0.4, 40)
  pr <- predict(mod, u, h = 1)
  ref <- cond_moments(qnorm(u), acf, 1)
  expect_equal(c(pr$mean, pr$sd), unname(ref), tolerance = 1e-8)
})

test_that("df, qf and dens are mutually consistent", {
  set.seed(8)
  mod <- armacop(ar = 0.5, ma = 0.4)
  u <- tssim(mod, 50)
  x <- c(0.05, 0.3, 0.5, 0.8, 0.97)
  for (h in c(1, 3)) {
    df <- predict(mod, u, h = h, type = "df", x = x)
    expect_equal(predict(mod, u, h = h, type = "qf", x = df), x, tolerance = 1e-8)
    expect_true(all(diff(df) > 0))
  }
  dens <- function(v) predict(mod, u, h = 2, type = "dens", x = v)
  expect_equal(integrate(dens, 0, 1)$value, 1, tolerance = 1e-5)
  eps <- 1e-6
  num <- (predict(mod, u, h = 2, type = "df", x = x + eps) -
            predict(mod, u, h = 2, type = "df", x = x - eps)) / (2 * eps)
  expect_equal(dens(x), num, tolerance = 1e-5)
})

test_that("predict works for fitted objects and validates input", {
  set.seed(9)
  u <- tssim(armacop(ar = 0.5), 200)
  f <- tsfit(u, armacop(order = c(1, 0)))
  expect_equal(predict(f), predict(f@model, u))
  expect_equal(predict(f, h = 2:3)$h, 2:3)
  expect_equal(predict(f, u[1:100]), predict(f@model, u[1:100]))
  mod <- f@model
  expect_error(predict(mod, u, h = 0), "'h'")
  expect_error(predict(mod, u, h = 1.5), "'h'")
  expect_error(predict(mod, u, h = 1:2, type = "df", x = 0.5), "single horizon")
  expect_error(predict(mod, u, type = "df"), "'x'")
  expect_error(predict(mod, u, type = "df", x = 1), "'x'")
  expect_error(predict(mod, u, type = "foo"), "arg")
  expect_error(predict(armacop(order = c(1, 0)), u), "unspecified")
})

test_that("predict agrees with tscopula one-step predictions", {
  skip_if_no_oracle()
  set.seed(10)
  mod <- armacop(ar = 0.5, ma = 0.4)
  u <- tssim(mod, 60)
  old <- tscopula::armacopula(list(ar = 0.5, ma = 0.4))
  x <- c(0.1, 0.5, 0.9)
  for (type in c("df", "qf", "dens")) {
    expect_equal(predict(mod, u, type = type, x = x),
                 as.numeric(tscopula::predict(old, u, x, type = type)), tolerance = 1e-6)
  }
})
