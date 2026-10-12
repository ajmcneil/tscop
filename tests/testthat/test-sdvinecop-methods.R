mixed_pairs <- function() {
  list(pair("gumbel", 90, 1.8), pair("clayton", 180, 1.2), pair("t", 0, c(0.3, 5)),
       pair("frank", 0, -2), pair("gauss", 0, 0.2), pair("bb1", 270, c(0.3, 1.3)))
}
mixed_model <- function() {
  sdvinecop(armacop(ar = -0.5, ma = 0.3),
           pairs = list(pc("gumbel"), pc("clayton", pos = 180, neg = 270), pc("t", aux = 5),
                        pc("frank"), pc("gumbel", pos = 180)), tautol = 1e-3)
}

# ---- simulation and Rosenblatt are mutually inverse -----------------------

test_that("dv_rosenblatt inverts dv_simulate", {
  set.seed(1)
  for (pairs in list(mixed_pairs(),
                     lapply(c(0.5, -0.3, 0.2), function(r) pair("gauss", 0, r)),
                     list(pair("clayton", 0, 2), pair("gumbel", 90, 1.5)))) {
    for (n in c(3, 12, 40)) {
      w <- matrix(runif(n * 4, 0.02, 0.98), n, 4)
      u <- dv_simulate(pairs, w)
      expect_true(all(u > 0 & u < 1))
      for (r in 1:4) expect_equal(dv_rosenblatt(pairs, u[, r]), w[, r], tolerance = 1e-6)
    }
  }
})

test_that("a series shorter than the effective lag is handled", {
  set.seed(2)
  pairs <- mixed_pairs()
  w <- matrix(runif(4), 4, 1)
  u <- dv_simulate(pairs, w)
  expect_equal(dv_rosenblatt(pairs, as.numeric(u)), as.numeric(w), tolerance = 1e-8)
  expect_equal(dv_simulate(pairs, matrix(0.4, 1, 1))[1, 1], 0.4)
})

test_that("Rosenblatt residuals equal numerically integrated conditional distributions", {
  set.seed(3)
  pairs <- mixed_pairs()
  K <- length(pairs)
  u <- as.numeric(dv_simulate(pairs, matrix(runif(12), 12, 1)))
  r <- dv_rosenblatt(pairs, u)
  expect_equal(r[1], u[1])
  for (t in c(2, 3, 5, K + 1, K + 2, 12)) {
    k <- min(t - 1, K)
    past <- u[(t - k):(t - 1)]
    lpast <- dv_loglik(past, pairs)
    dens <- function(y) vapply(y, function(yy) exp(dv_loglik(c(past, yy), pairs) - lpast), numeric(1))
    cdf <- integrate(dens, 0, u[t], rel.tol = 1e-9)$value
    expect_equal(r[t], cdf, tolerance = 1e-5)
  }
})

test_that("the first values have the stationary distribution (no burn-in)", {
  set.seed(4)
  pairs <- lapply(arma_pacf(c(0.5, 0.3), numeric(0), 2), function(a) pair("gauss", 0, a))
  X <- qnorm(t(dv_simulate(pairs, matrix(runif(3 * 20000), 3, 20000))))
  expect_equal(unname(cov(X)), toeplitz(acf_ref(c(0.5, 0.3), numeric(0), 2)),
               tolerance = 0.05, ignore_attr = TRUE)
  expect_lt(max(abs(colMeans(X))), 0.05)
})

test_that("simulated series have uniform margins, iid uniform residuals and the right lag-1 tau", {
  set.seed(5)
  m <- mixed_model()
  pairs <- dv_resolve(m)
  U <- dv_simulate(pairs, matrix(runif(300 * 30), 300, 30))
  expect_gt(ks.test(as.numeric(U), "punif")$p.value, 0.001)
  R <- apply(U, 2, function(v) dv_rosenblatt(pairs, v))
  expect_gt(ks.test(as.numeric(R), "punif")$p.value, 0.001)
  Z <- qnorm(R)
  expect_lt(abs(mean(Z)), 0.05)
  expect_equal(var(as.numeric(Z)), 1, tolerance = 0.05)
  for (lag in 1:4) {
    cors <- vapply(1:30, function(r) cor(Z[-(1:lag), r], Z[1:(300 - lag), r]), numeric(1))
    expect_lt(abs(mean(cors)), 0.03)
  }
  tau1 <- mean(vapply(1:30, function(r) cor(U[-1, r], U[-300, r], method = "kendall"), numeric(1)))
  expect_equal(tau1, ktau(m, 1), tolerance = 0.05)
})

test_that("tssim returns copula data and is reproducible", {
  m <- mixed_model()
  set.seed(6)
  a <- tssim(m, 60)
  set.seed(6)
  b <- tssim(m, 60)
  expect_equal(a, b)
  expect_length(a, 60)
  expect_true(all(a > 0 & a < 1))
  expect_length(tssim(sdvinecop(armacop(ar = 0.5)), 1), 1)
  expect_error(tssim(m, 0), "positive")
  expect_error(tssim(sdvinecop(armacop(order = c(1, 0)), pairs = list(pc("gumbel"))), 10), "unspecified")
  expect_error(tssim(sdvinecop(armacop(ar = 0.3), pairs = list(pc("bb1", aux = 3))), 10), "infeasible")
})

test_that("a Gaussian d-vine simulates the ARMA autocorrelation", {
  set.seed(7)
  m <- sdvinecop(armacop(ar = 0.6, ma = 0.3), tautol = 1e-6)
  e <- arma_coef_expanded(m@kpacf)
  x <- qnorm(tssim(m, 1500))
  expect_equal(as.numeric(acf(x, lag.max = 3, plot = FALSE)$acf[-1]),
               acf_ref(e$ar, e$ma, 3)[-1], tolerance = 0.1)
})

# ---- residuals -----------------------------------------------------------

test_that("resid agrees with armacop for a finite-order Gaussian model", {
  set.seed(8)
  arma <- armacop(ar = c(0.5, 0.3))
  m <- sdvinecop(arma)
  u <- tssim(arma, 80)
  expect_equal(resid(m, u, scale = "normal"), resid(arma, u, scale = "normal"), tolerance = 1e-8)
  expect_equal(resid(m, u), resid(arma, u), tolerance = 1e-8)
  expect_equal(residuals(m, u), resid(m, u))
})

test_that("resid is close to armacop for an ARMA copula at a small tolerance", {
  set.seed(9)
  arma <- armacop(ar = 0.5, ma = 0.4)
  m <- sdvinecop(arma, tautol = 1e-10)
  u <- tssim(arma, 80)
  expect_equal(resid(m, u, scale = "normal"), resid(arma, u, scale = "normal"), tolerance = 1e-6)
})

test_that("resid handles scale, attributes and invalid data", {
  set.seed(10)
  m <- mixed_model()
  u <- tssim(m, 40)
  expect_equal(resid(m, u, scale = "normal"), qnorm(resid(m, u)))
  ut <- ts(u, start = c(2000, 1), frequency = 4)
  expect_s3_class(resid(m, ut), "ts")
  expect_equal(tsp(resid(m, ut)), tsp(ut))
  u2 <- u
  u2[5] <- NA
  expect_error(resid(m, u2), "complete")
  expect_error(resid(m, c(0, u[-1])), "boundaryadjust")
  expect_error(resid(sdvinecop(armacop(order = c(1, 0)), pairs = list(pc("gumbel"))), u), "unspecified")
})

# ---- prediction ----------------------------------------------------------

test_that("predictive df equals the Rosenblatt residual of an appended value", {
  set.seed(11)
  m <- mixed_model()
  u <- tssim(m, 30)
  x <- c(0.05, 0.3, 0.5, 0.8, 0.97)
  df <- predict(m, u, type = "df", x = x)
  for (i in seq_along(x)) {
    expect_equal(df[i], resid(m, c(u, x[i]))[31], tolerance = 1e-8)
  }
})

test_that("predictive df, qf and dens are mutually consistent", {
  set.seed(12)
  m <- mixed_model()
  u <- tssim(m, 30)
  x <- c(0.05, 0.3, 0.5, 0.8, 0.97)
  df <- predict(m, u, type = "df", x = x)
  expect_true(all(diff(df) > 0))
  expect_equal(predict(m, u, type = "qf", x = df), x, tolerance = 1e-6)
  dens <- function(v) predict(m, u, type = "dens", x = v)
  expect_equal(integrate(dens, 0, 1, rel.tol = 1e-8)$value, 1, tolerance = 1e-4)
  eps <- 1e-6
  num <- (predict(m, u, type = "df", x = x + eps) - predict(m, u, type = "df", x = x - eps)) / (2 * eps)
  expect_equal(dens(x), num, tolerance = 1e-4)
})

test_that("predictions use the last K values and agree with armacop for a finite-order Gaussian model", {
  set.seed(13)
  arma <- armacop(ar = c(0.5, 0.3))
  m <- sdvinecop(arma)
  u <- tssim(arma, 50)
  x <- c(0.1, 0.5, 0.9)
  for (type in c("df", "qf", "dens")) {
    expect_equal(predict(m, u, type = type, x = x),
                 predict(arma, u, type = type, x = x), tolerance = 1e-8)
    expect_equal(predict(m, u[49:50], type = type, x = x), predict(m, u, type = type, x = x))
  }
  # data shorter than K
  mm <- mixed_model()
  expect_length(predict(mm, u[1:2], type = "df", x = x), 3)
  expect_length(predict(mm, u[1], type = "dens", x = x), 3)
})

test_that("predict validates its arguments", {
  m <- mixed_model()
  u <- tssim(m, 20)
  expect_error(predict(m, u, h = 2, type = "df", x = 0.5), "h = 1")
  expect_error(predict(m, u, type = "df"), "'x'")
  expect_error(predict(m, u, type = "df", x = 1), "'x'")
  expect_error(predict(m, u, type = "foo", x = 0.5), "arg")
  expect_error(predict(m, c(u, NA), type = "df", x = 0.5), "complete")
})

test_that("predict agrees with tscopula for a Gaussian-tail model", {
  skip_if_no_oracle()
  set.seed(14)
  arma <- armacop(ar = 0.5, ma = 0.4)
  m <- sdvinecop(arma, tautol = 1e-8)
  u <- tssim(arma, 80)
  old <- tscopula::armacopula(list(ar = 0.5, ma = 0.4))
  x <- c(0.1, 0.5, 0.9)
  for (type in c("df", "qf", "dens")) {
    expect_equal(predict(m, u, type = type, x = x),
                 as.numeric(tscopula::predict(old, u, x, type = type)), tolerance = 1e-5)
  }
})

# ---- ktau ----------------------------------------------------------------

test_that("ktau returns the ARMA kpacf up to K and zero beyond", {
  m <- mixed_model()
  k <- ktau(m, 30)
  expect_length(k, 30)
  e <- arma_coef_expanded(m@kpacf)
  ref <- (2 / pi) * asin(arma_pacf(e$ar, e$ma, m@K))
  expect_equal(k[seq_len(m@K)], ref)
  if (m@K < 30) expect_equal(k[(m@K + 1):30], rep(0, 30 - m@K))
  expect_equal(ktau(m, 2), ref[1:2])
  expect_length(ktau(m), 20)
  expect_error(ktau(m, 0), "lagmax")
})

test_that("ktau in free mode uses the Kendall's tau of each copula", {
  b <- list(rvinecopulib::bicop_dist("gumbel", 90, 2), rvinecopulib::bicop_dist("clayton", 0, 1))
  m <- sdvinecop(pairs = b)
  expect_equal(ktau(m, 4), c(-0.5, 1 / 3, 0, 0))
})

test_that("ktau of a Gaussian d-vine equals ktau of the armacop", {
  arma <- armacop(ar = 0.7, ma = 0.3)
  m <- sdvinecop(arma, tautol = 1e-12)
  expect_equal(ktau(m, 10), ktau(arma, 10), tolerance = 1e-8)
})
