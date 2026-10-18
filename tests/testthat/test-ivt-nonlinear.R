# Pure-R reference of nonlinear (v2p) v-transforms in IVT pair copulas (docs/sdvinecop-nonlinear-vt-spec.md).

nl_pair <- function(base = "joe", tau = 0.15, vt = c(0.62, 0.45), kappa = c(1.4, 0.7)) tscop:::dv_ivtnl_pair(base, tau, vt, kappa)

test_that("the v2p transform equals udp::v2p, is UDP, and kappa = 1 is the linear transform", {
  x <- seq(0.001, 0.999, length.out = 41)
  for (pr in list(c(0.3, 0.6), c(0.62, 1.7), c(0.45, 3))) {
    expect_equal(tscop:::dv_v2p(x, pr[1], pr[2]), as.numeric(udp::udptrans(udp::v2p(pr[1], pr[2]), x)), tolerance = 1e-12)
    g <- (seq_len(20000) - 0.5) / 20000
    tg <- tscop:::dv_v2p(g, pr[1], pr[2])
    expect_equal(c(mean(tg), var(tg)), c(0.5, 1 / 12), tolerance = 2e-3)   # T(U) is uniform
    h <- 1e-6
    xi <- c(0.1, 0.35, 0.8)
    expect_equal(tscop:::dv_v2p_slope(xi, pr[1], pr[2]),
                 (tscop:::dv_v2p(xi + h, pr[1], pr[2]) - tscop:::dv_v2p(xi - h, pr[1], pr[2])) / (2 * h), tolerance = 1e-6)
  }
  expect_equal(tscop:::dv_v2p(x, 0.4, 1), tscop:::dv_vlin(x, 0.4), tolerance = 1e-12)
})

test_that("the nonlinear IVT density is a copula density: integral one and uniform margins", {
  n <- 500
  g <- (seq_len(n) - 0.5) / n
  for (b in c("joe", "ast", "clayton180")) {
    pr <- nl_pair(b, 0.15)
    dens <- matrix(tscop:::dv_dens(pr, rep(g, times = n), rep(g, each = n)), n)
    expect_equal(mean(dens), 1, tolerance = 3e-3)
    expect_equal(rowMeans(dens)[c(40, 250, 450)], rep(1, 3), tolerance = 3e-3)
    expect_equal(colMeans(dens)[c(40, 250, 450)], rep(1, 3), tolerance = 3e-3)
  }
})

test_that("h-functions are distribution functions whose derivatives are the density, and the inverse inverts them", {
  set.seed(3)
  for (b in c("joe", "ast", "clayton180")) for (kp in list(c(1.4, 0.7), c(1, 0.7), c(1.6, 1))) {
    pr <- nl_pair(b, 0.2, kappa = kp)
    a <- runif(6, 0.03, 0.97)
    v <- runif(6, 0.03, 0.97)
    e <- 1e-5
    expect_equal((tscop:::dv_h(pr, a, v + e, 1L) - tscop:::dv_h(pr, a, v - e, 1L)) / (2 * e), tscop:::dv_dens(pr, a, v), tolerance = 2e-4)
    expect_equal((tscop:::dv_h(pr, a + e, v, 2L) - tscop:::dv_h(pr, a - e, v, 2L)) / (2 * e), tscop:::dv_dens(pr, a, v), tolerance = 2e-4)
    grid <- c(1e-9, seq(0.05, 0.95, 0.1), 1 - 1e-9)
    hh <- tscop:::dv_h(pr, 0.3, grid, 1L)
    expect_true(all(diff(hh) > 0))
    expect_lt(hh[1], 1e-6)
    expect_gt(hh[length(hh)], 1 - 1e-6)
    p <- runif(8, 0.02, 0.98)
    aa <- runif(8, 0.03, 0.97)
    expect_equal(tscop:::dv_h(pr, aa, tscop:::dv_hinv(pr, aa, p), 1L), p, tolerance = 1e-9)
  }
})

test_that("kappa = 1 reproduces the linear IVT pair copula and D-vine likelihood", {
  set.seed(8)
  tau <- as.numeric(ktau(arfactcop(lambda = 0.9, w = 0.5), 6))
  vt <- c(0.62, 0.45)
  u <- tssim(sdvinecop(arfactcop(lambda = 0.9, w = 0.5), tail = pc("ivt", base = "joe"), vt = vt, maxlag = 6), 400)
  lin <- lapply(tau, function(t) tscop:::dv_ivt_pair("joe", t, vt))
  nl <- tscop:::dv_nl_pairs("joe", tau, vt, c(1, 1), L = 6)
  a <- runif(5, 0.05, 0.95)
  b <- runif(5, 0.05, 0.95)
  expect_equal(tscop:::dv_dens(nl[[1]], a, b), tscop:::dv_dens(lin[[1]], a, b), tolerance = 1e-10)
  for (cv in 1:2) expect_equal(tscop:::dv_h(nl[[1]], a, b, cv), tscop:::dv_h(lin[[1]], a, b, cv), tolerance = 1e-10)
  expect_equal(tscop:::dv_loglik_r(u, nl), tscop:::dv_loglik(u, lin), tolerance = 1e-8)
  expect_null(tscop:::dv_encode(nl))   # nonlinear levels are not supported by the compiled kernel yet
  expect_equal(tscop:::dv_loglik(u, nl), tscop:::dv_loglik_r(u, nl))
})

test_that("the u-scale quadrature is converged and agrees with an independent probability-scale evaluation", {
  set.seed(9)
  pr <- nl_pair("joe", 0.15, c(0.62, 0.45), c(1.4, 0.7))
  a <- c(runif(20, 0.02, 0.98), 1e-3, 0.999)
  b <- c(runif(20, 0.02, 0.98), 0.999, 1e-3)
  ref <- function(cv, ...) {
    if (cv == 1L) tscop:::dv_nl_cdf(pr$base, pr$delta[2L], pr$kappa[2L], pr$delta[1L], pr$kappa[1L], b, a, FALSE, ...)
    else tscop:::dv_nl_cdf(pr$base, pr$delta[1L], pr$kappa[1L], pr$delta[2L], pr$kappa[2L], a, b, TRUE, ...)
  }
  for (cv in 1:2) {
    base24 <- ref(cv)
    expect_equal(base24, ref(cv, nodes = 96L, grade = 4), tolerance = 1e-4)
    expect_equal(base24[1:20], ref(cv, nodes = 96L, grade = 4)[1:20], tolerance = 2e-6)
  }
  # probability-scale form: h(x | g) = G(p0) + 1{x > d} p0, G(p0) = int_{p0}^1 Delta_x(Hinv(p | t_g)) dp (64-node rule)
  vtx <- udp::v2p(pr$delta[2L], pr$kappa[2L])
  tg <- tscop:::dv_v2p(a, pr$delta[1L], pr$kappa[1L])
  xs <- b[1:8]
  gl <- tscop:::gl_nodes(64)
  prob_scale <- vapply(1:8, function(i) {
    p0 <- rvinecopulib::hbicop(cbind(tscop:::dv_v2p(xs[i], pr$delta[2L], pr$kappa[2L]), tg[i]), 2L, "joe", 0, pr$base$parameters)
    P <- p0 + (1 - p0) * gl$s
    Q <- rvinecopulib::hbicop(cbind(P, rep(tg[i], 64)), 2L, "joe", 0, pr$base$parameters, inverse = TRUE)
    Dl <- attr(udp::udpinverse(vtx, pmin(pmax(as.numeric(Q), 1e-9), 1 - 1e-9), prob = TRUE), "prob")[, 1L]
    (1 - p0) * sum(gl$w * Dl) + (xs[i] > pr$delta[2L]) * p0
  }, numeric(1))
  expect_equal(ref(1L)[1:8], as.numeric(prob_scale), tolerance = 2e-4)
})

test_that("simulation inverts the Rosenblatt transform, and simulated residuals are uniform and independent", {
  set.seed(11)
  tau <- as.numeric(ktau(arfactcop(lambda = 0.9, w = 0.5), 8))
  pairs <- tscop:::dv_nl_pairs("joe", tau, c(0.62, 0.45), c(1.5, 0.75), L = 3)
  w <- matrix(runif(80), 80, 1)
  u <- tscop:::dv_simulate_r(pairs, w)
  expect_equal(tscop:::dv_rosenblatt_r(pairs, as.numeric(u)), as.numeric(w), tolerance = 1e-8)
  u2 <- as.numeric(tscop:::dv_simulate_r(pairs, matrix(runif(1500), 1500, 1)))
  r <- tscop:::dv_rosenblatt_r(pairs, u2)
  expect_gt(stats::ks.test(r[-(1:8)], "punif")$p.value, 0.001)
  z <- qnorm(r[-(1:8)])
  expect_lt(abs(cor(z[-1], z[-length(z)])), 0.1)
  expect_lt(abs(cor(abs(z[-1]), abs(z[-length(z)]))), 0.1)   # no remaining volatility clustering
})

test_that("the number of nonlinear lags follows the kpacf", {
  tau <- c(0.1, 0.05, 0.025, 0.015, 0.001)
  expect_equal(tscop:::dv_nl_lags(tau, 0.02), 3L)
  expect_equal(tscop:::dv_nl_lags(tau, 0.001), 4L)
  expect_equal(tscop:::dv_nl_lags(c(0.01, 0.005), 0.02), 1L)   # at least one nonlinear lag
  expect_equal(tscop:::dv_nl_lags(tau, 0), 5L)                 # capped at K
  pairs <- tscop:::dv_nl_pairs("ast", tau, c(0.6, 0.5), c(1.3, 0.8))
  expect_equal(vapply(pairs, function(p) p$family, ""), c("ivtnl", "ivtnl", "ivtnl", "ivt", "ivt"))
  expect_equal(tscop:::dv_nl_pairs("ast", tau, c(0.6, 0.5), c(1.3, 0.8), L = 2)[[3]]$family, "ivt")
})
