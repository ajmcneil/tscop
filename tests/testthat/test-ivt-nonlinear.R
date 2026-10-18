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

test_that("the compiled kernel evaluates nonlinear IVT levels like the R reference", {
  set.seed(21)
  tau <- as.numeric(ktau(arfactcop(lambda = 0.9, w = 0.5), 12))
  for (b in c("joe", "gumbel", "clayton180", "frank")) for (kp in list(c(1.3, 0.8), c(1, 0.8), c(0.7, 1), c(1, 1))) {
    pairs <- tscop:::dv_nl_pairs(b, tau, c(0.62, 0.45), kp, L = 5)
    u <- as.numeric(tscop:::dv_simulate_r(pairs, matrix(runif(500), 500, 1)))
    e <- tscop:::dv_encode(pairs, ast = TRUE, nl = TRUE)
    expect_false(is.null(e))
    expect_equal(tscop:::dv_loglik_cpp(u, e$fam, e$rot, e$par, e$ivt, e$bfam, e$brot, e$bpar, e$d1, e$d2, e$nl, e$k1, e$k2),
                 tscop:::dv_loglik_r(u, pairs), tolerance = 1e-8)
  }
  # nonlinear levels followed by a plain Gaussian run, a linear IVT run and plain levels (scale changes between runs)
  pairs <- c(tscop:::dv_nl_pairs("joe", tau[1:3], c(0.62, 0.45), c(1.3, 0.8), L = 3),
             list(list(family = "gauss", rotation = 0, parameters = 0.05)),
             list(tscop:::dv_ivt_pair("joe", tau[5], c(0.62, 0.45))),
             list(list(family = "gumbel", rotation = 0, parameters = 1.05)),
             list(tscop:::dv_ivtnl_pair("joe", tau[7], c(0.62, 0.45), c(1.3, 0.8))))
  u <- runif(400)
  expect_equal(tscop:::dv_loglik(u, pairs), tscop:::dv_loglik_r(u, pairs), tolerance = 1e-8)
})

test_that("the compiled nonlinear IVT gradient equals finite differences, including the kappa directions", {
  set.seed(22)
  tau <- as.numeric(ktau(arfactcop(lambda = 0.95, w = 0.25), 16))
  for (b in c("joe", "clayton180", "frank")) for (kp in list(c(1.3, 0.8), c(1, 1))) {
    pairs <- tscop:::dv_nl_pairs(b, tau, c(0.62, 0.45), kp, L = 6)
    u <- as.numeric(tscop:::dv_simulate_r(tscop:::dv_nl_pairs(b, tau, c(0.62, 0.45), c(1.3, 0.8), L = 6), matrix(runif(400), 400, 1)))
    e <- tscop:::dv_encode(pairs, ast = TRUE, nl = TRUE)
    K <- length(e$fam)
    dbpar <- matrix(0, K, 5)
    dbpar[, 5] <- e$bpar                       # fifth direction: relative change of the base parameters
    g <- tscop:::dv_loglik_grad_cpp(u, matrix(0, 0, 5), e$fam, e$rot, e$par, matrix(0, K, 5), e$ivt, e$bfam, e$brot, e$bpar, dbpar,
                                    e$d1, e$d2, c(1, 0, 0, 0, 0), c(0, 1, 0, 0, 0), e$nl, e$k1, e$k2, c(0, 0, 1, 0, 0), c(0, 0, 0, 1, 0))
    f <- function(x) tscop:::dv_loglik_cpp(u, e$fam, e$rot, e$par, e$ivt, e$bfam, e$brot, e$bpar * (1 + x[5]), x[1], x[2], e$nl, x[3], x[4])
    x0 <- c(e$d1, e$d2, e$k1, e$k2, 0)
    h <- 1e-5
    fd <- vapply(1:5, function(i) { xp <- x0; xp[i] <- xp[i] + h; xm <- x0; xm[i] <- xm[i] - h; (f(xp) - f(xm)) / (2 * h) }, 1)
    expect_lt(max(abs(g$gradient - fd) / pmax(abs(fd), 1)), 2e-4)
    expect_equal(g$value, f(x0), tolerance = 1e-6)
  }
})

test_that("nonlinear levels are encoded only when the caller handles them and not for the Gaussian base", {
  tau <- as.numeric(ktau(arfactcop(lambda = 0.9, w = 0.5), 5))
  p <- tscop:::dv_nl_pairs("joe", tau, c(0.6, 0.5), c(1.3, 0.8), L = 3)
  expect_null(tscop:::dv_encode(p, ast = TRUE))                 # callers that do not handle nonlinear levels (ktauhat) use the R path
  expect_false(is.null(tscop:::dv_encode(p, ast = TRUE, nl = TRUE)))
  expect_false(is.null(tscop:::dv_encode(tscop:::dv_nl_pairs("ast", tau, c(0.6, 0.5), c(1.3, 0.8), L = 3), ast = TRUE, nl = TRUE)))
  pg <- tscop:::dv_nl_pairs("gauss", tau, c(0.6, 0.5), c(1.3, 0.8), L = 3)
  expect_null(tscop:::dv_encode(pg, ast = TRUE, nl = TRUE))     # Gaussian base: R path
  u <- runif(300)
  expect_true(is.finite(tscop:::dv_loglik(u, pg)))
  expect_equal(tscop:::dv_loglik(u, pg), tscop:::dv_loglik_r(u, pg))
})

test_that("the ast quantile table is accurate", {
  v <- c(1e-10, 1e-6, 1e-3, 0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99, 0.999, 1 - 1e-5, 1 - 1e-8, 1 - 1e-10)
  for (nu in c(1.2, 3, 20, 200)) {
    tb <- tscop:::nl_ast_table_cpp(v, nu)
    z <- qt((1 - v) / 2, nu, lower.tail = FALSE)
    expect_lt(max(abs(tb[, 1] - z) / pmax(z, 1e-3)), 1e-4)
    expect_lt(max(abs(tb[v > 1e-8, 2] * 2 * dt(z[v > 1e-8], nu) - 1)), 2e-4)   # dz/dV = 1 / (2 f_t(z))
    h <- 1e-4 * nu
    dn <- (qt((1 - v) / 2, nu + h, lower.tail = FALSE) - qt((1 - v) / 2, nu - h, lower.tail = FALSE)) / (2 * h)
    expect_lt(max(abs(tb[v > 0.1, 3] - dn[v > 0.1]) / abs(dn[v > 0.1])), 1e-4)   # the nu-derivative matters where z is not small
    expect_lt(max(abs(tb[, 3] - dn)) / max(abs(dn)), 1e-5)
  }
})

test_that("ast bases of nonlinear levels: compiled kernels agree with the R reference on model data, and the gradient with finite differences", {
  set.seed(41)
  tau <- as.numeric(ktau(arfactcop(lambda = 0.95, w = 0.25), 14))
  for (kp in list(c(1.3, 0.8), c(1, 0.8))) {
    pairs <- tscop:::dv_nl_pairs("ast", tau, c(0.62, 0.45), kp, L = 6)
    u <- as.numeric(tscop:::dv_simulate(pairs, matrix(runif(500), 500, 1)))
    expect_equal(tscop:::dv_loglik(u, pairs), tscop:::dv_loglik_r(u, pairs), tolerance = 1e-6)
    expect_equal(tscop:::dv_rosenblatt(pairs, u), tscop:::dv_rosenblatt_r(pairs, u), tolerance = 1e-4)
    x <- c(0.02, 0.2, 0.5, 0.8, 0.98)
    for (ty in c("df", "qf", "dens")) expect_equal(tscop:::dv_predictive(pairs, u, x, ty), tscop:::dv_predictive_r(pairs, u, x, ty), tolerance = 1e-5)
    s_r <- tscop:::dv_simulate_r(pairs, matrix(runif(100), 100, 1))
    expect_equal(tscop:::dv_simulate(pairs, matrix(runif(100), 100, 1)) > 0, matrix(TRUE, 100, 1))
    e <- tscop:::dv_encode(pairs, ast = TRUE, nl = TRUE)
    K <- length(e$fam)
    dbpar <- matrix(0, K, 5)
    dbpar[, 5] <- e$bpar
    g <- tscop:::dv_loglik_grad_cpp(u, matrix(0, 0, 5), e$fam, e$rot, e$par, matrix(0, K, 5), e$ivt, e$bfam, e$brot, e$bpar, dbpar,
                                    e$d1, e$d2, c(1, 0, 0, 0, 0), c(0, 1, 0, 0, 0), e$nl, e$k1, e$k2, c(0, 0, 1, 0, 0), c(0, 0, 0, 1, 0))
    f <- function(z) tscop:::dv_loglik_cpp(u, e$fam, e$rot, e$par, e$ivt, e$bfam, e$brot, e$bpar * (1 + z[5]), z[1], z[2], e$nl, z[3], z[4])
    x0 <- c(e$d1, e$d2, e$k1, e$k2, 0)
    h <- 1e-5
    fd <- vapply(1:5, function(i) { xp <- x0; xp[i] <- xp[i] + h; xm <- x0; xm[i] <- xm[i] - h; (f(xp) - f(xm)) / (2 * h) }, 1)
    expect_lt(max(abs(g$gradient - fd) / pmax(abs(fd), 1)), 1e-4)
  }
})
