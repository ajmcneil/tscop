# The compiled likelihood kernel against the R reference implementation.

cmp_cpp <- function(pairs, u, tol = 1e-9) {
  r <- dv_loglik_r(u, pairs)
  cpp <- dv_loglik(u, pairs)
  expect_false(is.null(dv_encode(pairs)))
  expect_equal(cpp, r, tolerance = tol)
}

test_that("the C++ kernel equals the R likelihood for plain families", {
  set.seed(1)
  kp <- arfactcop(lambda = 0.95, w = 0.5)
  for (f in c("joe", "gumbel", "clayton", "frank", "gauss")) {
    m <- sdvinecop(kp, tail = pc(f), maxlag = 30, tautol = 1e-6)
    m@K <- 30L
    cmp_cpp(dv_resolve(m), tssim(m, 800))
  }
})

test_that("the C++ kernel equals the R likelihood for rotations and negative dependence", {
  set.seed(2)
  arm <- armacop(ar = c(0.5, -0.3))
  for (rot in list(c(0, 90), c(180, 270))) {
    m <- sdvinecop(arm, pairs = list(pc("gumbel", pos = rot[1], neg = rot[2]), pc("clayton", pos = rot[1], neg = rot[2]),
                                     pc("joe", pos = rot[1], neg = rot[2]), pc("frank")), tail = "gauss", maxlag = 15)
    m@K <- 15L
    pr <- dv_resolve(m)
    # the first three lags only: the later lags have Kendall's tau near zero, where
    # rvinecopulib's Frank density is itself inaccurate (see the next test)
    cmp_cpp(pr[1:3], tssim(m, 800), tol = 1e-7)
  }
  m <- sdvinecop(armacop(ar = -0.6, ma = 0.2), pairs = list(pc("gumbel"), pc("joe"), pc("clayton"), pc("frank")),
                 tail = pc("gumbel"), maxlag = 12)
  m@K <- 12L
  cmp_cpp(dv_resolve(m), tssim(m, 800))
})

test_that("the C++ Frank kernel is accurate near independence where rvinecopulib is not", {
  set.seed(3)
  u <- runif(2000)
  pr <- list(family = "frank", rotation = 0, parameters = 1e-10)
  # the true log-likelihood is of order theta; rvinecopulib's density has errors of about 4e-6 per point
  expect_lt(abs(dv_loglik(u, list(pr))), 1e-8)
  expect_gt(abs(dv_loglik_r(u, list(pr))), 1e-5)
  th <- 0.05
  expect_equal(dv_loglik(u, list(modifyList(pr, list(parameters = th)))),
               dv_loglik_r(u, list(modifyList(pr, list(parameters = th)))), tolerance = 1e-9)
})

test_that("the C++ kernel equals the R likelihood for IVT pairs with every base", {
  set.seed(4)
  kp <- arfactcop(lambda = 0.95, w = 0.5)
  for (b in c("joe", "gumbel", "clayton180", "frank", "gauss")) {
    for (d in list(c(0.5, 0.5), c(0.3, 0.7))) {
      m <- sdvinecop(kp, tail = pc("ivt", base = b), vt = d, maxlag = 25, tautol = 1e-6)
      m@K <- 25L
      cmp_cpp(dv_resolve(m), tssim(m, 600))
    }
  }
  mm <- sdvinecop(kp, pairs = c(rep(pc("ivt", base = "joe"), 3), list(pc("gumbel")), rep(pc("ivt", base = "frank"), 2)),
                  tail = pc("ivt", base = "gauss"), vt = c(0.4, 0.6), maxlag = 20)
  mm@K <- 20L
  cmp_cpp(dv_resolve(mm), tssim(mm, 600))
  m2 <- sdvinecop(kp, pairs = c(list(pc("gauss")), rep(pc("ivt", base = "joe"), 2), rep(pc("gauss"), 3)),
                  tail = pc("gumbel"), vt = c(0.45, 0.55), maxlag = 15)
  m2@K <- 15L
  cmp_cpp(dv_resolve(m2), runif(500))
})

test_that("unsupported pair copulas use the R implementation and the option switches the kernel off", {
  kp <- arfactcop(lambda = 0.9, w = 0.5)
  m <- unsupported_model()
  pr <- dv_resolve(m)
  expect_null(dv_encode(pr, ast = TRUE))
  u <- tssim(m, 300)
  expect_equal(dv_loglik(u, pr), dv_loglik_r(u, pr))
  ma <- sdvinecop(kp, tail = pc("ivt", base = "ast"), maxlag = 5)
  ma@K <- 5L
  expect_null(dv_encode(dv_resolve(ma)))
  mj <- sdvinecop(kp, tail = pc("joe"), maxlag = 8)
  mj@K <- 8L
  uj <- tssim(mj, 300)
  old <- options(tscop.cpp = FALSE)
  on.exit(options(old))
  expect_equal(dv_loglik(uj, dv_resolve(mj)), dv_loglik_r(uj, dv_resolve(mj)))
})

test_that("models fitted with the kernel agree with fits by the R implementation", {
  skip_on_cran()
  set.seed(5)
  truth <- sdvinecop(arfactcop(lambda = 0.9, w = 0.5), pairs = list(pc("gumbel")), tail = "gauss", maxlag = 15)
  u <- tssim(truth, 600)
  spec <- sdvinecop(arfactcop(order = 1), pairs = list(pc("gumbel")), tail = "gauss", maxlag = 15)
  f1 <- tsfit(u, spec, pseudo = FALSE)
  old <- options(tscop.cpp = FALSE)
  on.exit(options(old))
  f2 <- tsfit(u, spec, pseudo = FALSE)
  expect_equal(unname(coef(f1)), unname(coef(f2)), tolerance = 1e-4)
  expect_equal(as.numeric(logLik(f1)), as.numeric(logLik(f2)), tolerance = 1e-7)
})

# ---- exact gradients -------------------------------------------------------

test_that("the kernel gradient with respect to the level parameters and fulcrums matches finite differences", {
  set.seed(6)
  K <- 8L
  kp <- arfactcop(lambda = 0.95, w = 0.5)
  for (b in c("gumbel", "clayton180", "gauss")) {
    m <- sdvinecop(kp, tail = pc("ivt", base = b), vt = c(0.4, 0.6), maxlag = K, tautol = 1e-8)
    m@K <- K
    u <- tssim(m, 500)
    e <- dv_encode(dv_resolve(m))
    P <- K + 2L
    dpar <- dbp <- matrix(0, K, P)
    for (j in seq_len(K)) dbp[j, j] <- 1
    dd1 <- dd2 <- numeric(P)
    dd1[K + 1L] <- 1
    dd2[K + 2L] <- 1
    g <- dv_loglik_grad_cpp(u, matrix(0, 0, ncol(dpar)), e$fam, e$rot, e$par, dpar, e$ivt, e$bfam, e$brot, e$bpar, dbp, e$d1, e$d2, dd1, dd2)
    f <- function(bpar, d1, d2) dv_loglik_cpp(u, e$fam, e$rot, e$par, e$ivt, e$bfam, e$brot, bpar, d1, d2)
    expect_equal(g$value, f(e$bpar, e$d1, e$d2), tolerance = 1e-12)
    h <- 1e-7
    fd <- numeric(P)
    for (j in seq_len(K)) {
      bp <- e$bpar; bp[j] <- bp[j] + h; a <- f(bp, e$d1, e$d2)
      bp <- e$bpar; bp[j] <- bp[j] - h; fd[j] <- (a - f(bp, e$d1, e$d2)) / (2 * h)
    }
    fd[K + 1L] <- (f(e$bpar, e$d1 + h, e$d2) - f(e$bpar, e$d1 - h, e$d2)) / (2 * h)
    fd[K + 2L] <- (f(e$bpar, e$d1, e$d2 + h) - f(e$bpar, e$d1, e$d2 - h)) / (2 * h)
    expect_lt(max(abs(g$gradient - fd) / pmax(abs(fd), 1)), 2e-4)
  }
  # plain families, with rotations (the first three lags: a Frank pair at tau near zero has a
  # parameter of order 1e-10, which a finite-difference step of 1e-7 would cross)
  m <- sdvinecop(armacop(ar = c(0.5, -0.3)), pairs = list(pc("gumbel"), pc("clayton"), pc("joe")),
                 tail = pc("gumbel"), maxlag = K, tautol = 1e-3)
  m@K <- 3L
  K <- 3L
  u <- tssim(m, 500)
  e <- dv_encode(dv_resolve(m))
  P <- K
  dpar <- diag(K)
  dbp <- matrix(0, K, P)
  z <- numeric(P)
  g <- dv_loglik_grad_cpp(u, matrix(0, 0, ncol(dpar)), e$fam, e$rot, e$par, dpar, e$ivt, e$bfam, e$brot, e$bpar, dbp, e$d1, e$d2, z, z)
  f <- function(p) dv_loglik_cpp(u, e$fam, e$rot, p, e$ivt, e$bfam, e$brot, e$bpar, e$d1, e$d2)
  h <- 1e-7
  fd <- vapply(seq_len(K), function(j) { p <- e$par; p[j] <- p[j] + h; a <- f(p); p[j] <- p[j] - 2 * h; (a - f(p)) / (2 * h) }, 1)
  expect_lt(max(abs(g$gradient - fd) / pmax(abs(fd), 1)), 2e-4)
})

test_that("dv_negll_grad equals finite differences of dv_negll on the free scale", {
  set.seed(7)
  n <- 700
  K <- 20L
  cases <- list(
    sdvinecop(arfactcop(lambda = 0.95, w = 0.5), tail = pc("ivt", base = "gumbel"), vt = c(0.4, 0.6), maxlag = K, tautol = 1e-6),
    sdvinecop(arfactcop(lambda = 0.95, w = 0.5), tail = pc("ivt", base = "clayton180"), vt = c(0.4, 0.4),
              vtpar = "common", maxlag = K, tautol = 1e-6),
    sdvinecop(arfactcop(lambda = c(0.98, 0.6), w = c(0.3, 0.3)), tail = pc("ivt", base = "gauss"), vt = c(0.5, 0.5),
              maxlag = K, tautol = 1e-6),
    sdvinecop(armacop(ar = 0.7, ma = 0.2), pairs = list(pc("gumbel")), tail = pc("clayton"), maxlag = K, tautol = 1e-6),
    sdvinecop(arfactcop(lambda = 0.95, w = 0.5), tail = pc("joe"), maxlag = K, tautol = 1e-6),
    sdvinecop(arfactcop(lambda = 0.95, w = 0.5), tail = pc("ivt", base = "joe"), vt = c(0.4, 0.6), maxlag = K, tautol = 1e-6),
    sdvinecop(arfactcop(lambda = 0.95, w = 0.5), tail = pc("ivt", base = "frank"), vt = c(0.45, 0.55), maxlag = K, tautol = 1e-6)
  )
  for (m in cases) {
    m@K <- K
    u <- tssim(m, n)
    lay <- dv_layout(m)
    z0 <- dv_start(m, NULL, u, list())$z0
    z <- z0 + rnorm(length(z0), 0, 0.05)
    g <- dv_negll_grad(z, lay, K, u)
    expect_false(is.null(g))
    fd <- vapply(seq_along(z), function(k) {
      h <- 1e-6
      zp <- z; zp[k] <- zp[k] + h
      zm <- z; zm[k] <- zm[k] - h
      (dv_negll(zp, lay, K, u) - dv_negll(zm, lay, K, u)) / (2 * h)
    }, 1)
    expect_lt(max(abs(g - fd) / pmax(abs(fd), 1)), 1e-4)
  }
})

test_that("dv_negll_grad is exact for a pure AR kpacf, whose deep lags have tau of order 1e-17", {
  set.seed(7)
  K <- 6L
  m <- sdvinecop(armacop(ar = c(0.6, -0.3)), pairs = list(pc("joe"), pc("frank"), pc("joe")), tail = pc("frank"),
                 maxlag = K, tautol = 1e-6)
  m@K <- K
  u <- tssim(m, 700)
  lay <- dv_layout(m)
  z <- dv_start(m, NULL, u, list())$z0 + c(0.03, -0.02)
  g <- dv_negll_grad(z, lay, K, u)
  expect_false(is.null(g))
  # the lags beyond the AR order have tau of either sign, flipping the Frank parameter at 1e-10: a
  # finite-difference step of 1e-4 is above that noise
  fd <- vapply(seq_along(z), function(k) {
    h <- 1e-4
    zp <- z; zp[k] <- zp[k] + h
    zm <- z; zm[k] <- zm[k] - h
    (dv_negll(zp, lay, K, u) - dv_negll(zm, lay, K, u)) / (2 * h)
  }, 1)
  expect_lt(max(abs(g - fd) / pmax(abs(fd), 1)), 1e-5)
})

test_that("dv_negll_grad declines unsupported models and dv_negll_gr then falls back to differences", {
  set.seed(8)
  m <- unsupported_model()
  u <- tssim(m, 300)
  lay <- dv_layout(m)
  z <- dv_start(m, NULL, u, list())$z0
  expect_null(dv_negll_grad(z, lay, 6L, u))
  gr <- dv_negll_gr(z, lay, 6L, u)
  expect_length(gr, length(z))
  expect_true(all(is.finite(gr)))
})

# ---- gradients for v-transform copulas -----------------------------------------

test_that("vt_negll_grad equals finite differences for v-transform copulas with a d-vine base", {
  set.seed(9)
  kp <- arfactcop(lambda = 0.9, w = 0.5)
  b <- function(f) sdvinecop(kp, pairs = list(pc(f)), tail = pc(f), maxlag = 10, tautol = 1e-6)
  cases <- list(
    vtscop(b("joe"), udp::v2p(0.45, 0.8)),
    vtscop(b("gumbel"), udp::v2b(0.45, 0.8)),
    vtscop(b("clayton"), udp::vlinear(0.45)),
    vtscop(b("frank"), udp::v3p(0.45, 0.8, 1.2)),
    vtscop(b("joe"), udp::v2p(0.45, 0.8), vtfix = "delta")
  )
  for (m in cases) {
    u <- tssim(m, 500)
    L <- vt_layout(m)
    K <- 10L
    st <- vt_start(m, NULL, u, list())
    z <- st$z0 + rnorm(length(st$z0), 0, 0.03)
    g <- vt_negll_grad(z, L, K, u)
    expect_false(is.null(g))
    fd <- vapply(seq_along(z), function(k) {
      h <- 1e-6
      zp <- z; zp[k] <- zp[k] + h
      zm <- z; zm[k] <- zm[k] - h
      (vt_negll(zp, L, K, u) - vt_negll(zm, L, K, u)) / (2 * h)
    }, 1)
    expect_lt(max(abs(g - fd) / pmax(abs(fd), 1)), 1e-4)
  }
})

test_that("vt_negll_grad declines a missing-data series and an NA effective lag for a d-vine base", {
  set.seed(10)
  m <- vtscop(arfactcop(lambda = 0.9, w = 0.5), udp::v2p(0.45, 0.8))
  u <- tssim(m, 300)
  L <- vt_layout(m)
  z <- vt_start(m, NULL, u, list())$z0
  expect_false(is.null(vt_negll_grad(z, L, NA_integer_, u)))
  u2 <- u
  u2[5] <- NA
  expect_null(vt_negll_grad(z, L, NA_integer_, u2))
  expect_true(all(is.finite(vt_negll_gr(z, L, NA_integer_, u))))
  # a d-vine base needs the effective lag
  md <- vtscop(sdvinecop(arfactcop(lambda = 0.9, w = 0.5), pairs = list(pc("joe")), tail = pc("joe"), maxlag = 6), udp::v2p(0.45, 0.8))
  Ld <- vt_layout(md)
  zd <- vt_start(md, NULL, tssim(md, 200), list())$z0
  expect_null(vt_negll_grad(zd, Ld, NA_integer_, tssim(md, 200)))
})

test_that("vtscop fits with the exact gradient agree with fits by finite differences", {
  skip_on_cran()
  set.seed(11)
  truth <- vtscop(sdvinecop(arfactcop(lambda = 0.9, w = 0.5), pairs = list(pc("joe")), tail = pc("joe"), maxlag = 12, tautol = 1e-4),
                  udp::v2p(0.45, 0.8))
  u <- tssim(truth, 800)
  spec <- vtscop(sdvinecop(arfactcop(order = 1), pairs = list(pc("joe")), tail = pc("joe"), maxlag = 12, tautol = 1e-4),
                 udp::v2p(0.5, 1))
  f1 <- tsfit(u, spec, pseudo = FALSE)
  old <- options(tscop.cppgrad = FALSE)
  on.exit(options(old))
  f2 <- tsfit(u, spec, pseudo = FALSE)
  expect_equal(unname(coef(f1)), unname(coef(f2)), tolerance = 1e-3)
  expect_equal(as.numeric(logLik(f1)), as.numeric(logLik(f2)), tolerance = 1e-5)
})

# ---- the arfactcop Kalman filter ---------------------------------------------

test_that("the compiled Kalman filter equals the R filter, with and without missing values", {
  set.seed(14)
  for (cfg in list(list(l = 0.9, w = 0.6), list(l = c(0.98, 0.7), w = c(0.3, 0.4)), list(l = c(0.98, 0.8, 0.5), w = c(0.2, 0.2, 0.3)))) {
    g <- arfactcop(lambda = cfg$l, w = cfg$w)
    x <- stats::qnorm(tssim(g, 600))
    p <- list(lambda = cfg$l, w = cfg$w, w0 = 1 - sum(cfg$w))
    expect_equal(afc_loglik_cpp(x, p$lambda, p$w, p$w0, 1e-12), afc_filter(x, p$lambda, p$w, p$w0)$loglik, tolerance = 1e-12)
    x[c(1, 7, 8, 300, 600)] <- NA
    expect_equal(afc_loglik_cpp(x, p$lambda, p$w, p$w0, 1e-12), afc_filter(x, p$lambda, p$w, p$w0)$loglik, tolerance = 1e-12)
  }
})

test_that("afc_negll_grad equals finite differences, for parameters and for data directions", {
  set.seed(15)
  fd <- function(f, z, h = 1e-6) vapply(seq_along(z), function(k) {
    zp <- z; zp[k] <- zp[k] + h
    zm <- z; zm[k] <- zm[k] - h
    (f(zp) - f(zm)) / (2 * h)
  }, 1)
  for (cfg in list(list(order = 1, noise = TRUE), list(order = 2, noise = TRUE), list(order = 3, noise = TRUE),
                   list(order = 2, noise = FALSE))) {
    w <- if (cfg$noise) c(0.3, 0.2, 0.2)[seq_len(cfg$order)] else c(0.4, 0.6)
    g <- arfactcop(lambda = c(0.97, 0.8, 0.5)[seq_len(cfg$order)], w = w, noise = cfg$noise)
    x <- stats::qnorm(tssim(g, 800))
    x[c(3, 200, 201)] <- NA
    z <- gc_free_ops(g)$to_free(g) + rnorm(afc_nfree(cfg$order, cfg$noise), 0, 0.05)
    ex <- afc_negll_grad(z, cfg$order, cfg$noise, x)
    nu <- fd(function(zz) afc_negll(zz, cfg$order, cfg$noise, x), z)
    expect_lt(max(abs(ex - nu) / pmax(abs(nu), 1)), 1e-5)
  }
  g <- arfactcop(lambda = c(0.97, 0.7), w = c(0.3, 0.4))
  x0 <- stats::qnorm(tssim(g, 500))
  d1 <- sin(seq_along(x0) / 7)
  d2 <- cos(seq_along(x0) / 11) / 2
  z <- gc_free_ops(g)$to_free(g) + rnorm(4, 0, 0.05)
  f <- function(par) afc_negll(par[3:6], 2, TRUE, x0 + par[1] * d1 + par[2] * d2)
  par0 <- c(0.1, -0.2, z)
  ex <- afc_negll_grad(z, 2, TRUE, x0 + 0.1 * d1 - 0.2 * d2, cbind(d1, d2, matrix(0, 500, 4)), before = 2L)
  expect_lt(max(abs(ex - fd(f, par0)) / pmax(abs(fd(f, par0)), 1)), 1e-5)
})

test_that("the exact joint gradient also covers arfactcop and vtscop-on-arfactcop copulas", {
  set.seed(16)
  kp <- arfactcop(lambda = 0.9, w = 0.5)
  cases <- list(
    list(tsmod(marg("st", c(mu = 0, sigma = 1, df = 6)), kp), tsmod("st", arfactcop(order = 1))),
    list(tsmod(marg("sst", c(mu = 0, sigma = 1, df = 6, gamma = 1.2)), arfactcop(lambda = c(0.95, 0.6), w = c(0.3, 0.3))),
         tsmod("sst", arfactcop(order = 2))),
    list(tsmod(marg("st", c(mu = 0, sigma = 1, df = 8)), vtscop(kp, udp::v2p(0.45, 0.8))),
         tsmod("st", vtscop(arfactcop(order = 1), udp::v2p(0.5, 1))))
  )
  for (cs in cases) {
    x <- tssim(cs[[1]], 500)
    f2 <- tsfit(x, cs[[2]])
    jt <- tm_joint_setup(f2@model@margin, f2@model@copula, x)
    # (a fit can sit on the bound of the free scale, so shrink the point inwards before perturbing)
    th <- 0.98 * jt$theta0 + rnorm(length(jt$theta0), 0, 0.01)
    g <- jt$gradient(th, x)
    expect_false(is.null(g))
    fd <- vapply(seq_along(th), function(k) {
      h <- 1e-6
      tp <- th; tp[k] <- tp[k] + h
      tm <- th; tm[k] <- tm[k] - h
      (jt$negll(tp, x) - jt$negll(tm, x)) / (2 * h)
    }, 1)
    expect_lt(max(abs(g - fd) / pmax(abs(fd), 1)), 1e-4)
  }
})

# ---- the ARMA Kalman filter ------------------------------------------------------

test_that("the compiled ARMA filter equals the R filter for plain and seasonal models, with missing values", {
  set.seed(18)
  cfgs <- list(list(ar = 0.6, ma = numeric(0)), list(ar = numeric(0), ma = 0.5), list(ar = c(0.5, -0.2), ma = 0.4),
               list(ar = 0.7, ma = c(0.3, 0.1, -0.2)), list(ar = numeric(0), ma = numeric(0)),
               arma_expand(ar = 0.5, ma = 0.3, sar = 0.6, sma = 0.4, period = 12))
  for (cf in cfgs) {
    x <- stats::rnorm(300)
    expect_equal(arma_loglik_cpp(x, cf$ar, cf$ma, 1e-12), arma_filter(x, cf$ar, cf$ma)$loglik, tolerance = 1e-11)
    x[c(1, 4, 5, 150, 300)] <- NA
    expect_equal(arma_loglik_cpp(x, cf$ar, cf$ma, 1e-12), arma_filter(x, cf$ar, cf$ma)$loglik, tolerance = 1e-11)
  }
})

test_that("armacop_negll_grad equals finite differences, for coefficients and for data directions", {
  set.seed(19)
  fd <- function(f, z, h = 1e-6) vapply(seq_along(z), function(k) {
    zp <- z; zp[k] <- zp[k] + h
    zm <- z; zm[k] <- zm[k] - h
    (f(zp) - f(zm)) / (2 * h)
  }, 1)
  for (ord in list(c(1, 0, 0, 0), c(0, 1, 0, 0), c(1, 1, 0, 0), c(2, 1, 0, 0), c(1, 2, 0, 0), c(1, 1, 1, 1))) {
    period <- if (ord[3] + ord[4] > 0) 12L else NA
    x <- stats::rnorm(400)
    x[c(5, 6, 200)] <- NA
    z <- stats::rnorm(sum(ord), 0, 0.5)
    g <- armacop_negll_grad(z, ord, period, x)
    f <- fd(function(zz) armacop_negll(zz, ord, period, x), z)
    expect_lt(max(abs(g - f) / pmax(abs(f), 1)), 1e-6)
  }
  x0 <- stats::rnorm(300)
  d1 <- sin(seq_along(x0) / 7)
  d2 <- cos(seq_along(x0) / 11) / 2
  ord <- c(2, 1, 0, 0)
  z <- stats::rnorm(3, 0, 0.5)
  f <- function(par) armacop_negll(par[3:5], ord, NA, x0 + par[1] * d1 + par[2] * d2)
  par0 <- c(0.1, -0.2, z)
  g <- armacop_negll_grad(z, ord, NA, x0 + 0.1 * d1 - 0.2 * d2, cbind(d1, d2, matrix(0, 300, 3)), before = 2L)
  expect_lt(max(abs(g - fd(f, par0)) / pmax(abs(fd(f, par0)), 1)), 1e-6)
})

test_that("the exact joint gradient covers an armacop copula, and fits agree with differences", {
  skip_on_cran()
  set.seed(17)
  x <- tssim(tsmod(marg("st", c(mu = 0, sigma = 1, df = 6)), armacop(ar = 0.6, ma = 0.3)), 500)
  spec <- tsmod("st", armacop(order = c(1, 1)))
  f <- tsfit(x, spec)
  jt <- tm_joint_setup(f@model@margin, f@model@copula, x)
  th <- jt$theta0 + rnorm(length(jt$theta0), 0, 0.02)
  g <- jt$gradient(th, x)
  expect_false(is.null(g))
  fd <- vapply(seq_along(th), function(k) {
    h <- 1e-6
    tp <- th; tp[k] <- tp[k] + h
    tm <- th; tm[k] <- tm[k] - h
    (jt$negll(tp, x) - jt$negll(tm, x)) / (2 * h)
  }, 1)
  expect_lt(max(abs(g - fd) / pmax(abs(fd), 1)), 1e-5)
  f1 <- tsfit(x, spec, method = "joint")
  old <- options(tscop.cppgrad = FALSE)
  on.exit(options(old))
  f2 <- tsfit(x, spec, method = "joint")
  expect_equal(unname(coef(f1)), unname(coef(f2)), tolerance = 1e-3)
  expect_equal(as.numeric(logLik(f1)), as.numeric(logLik(f2)), tolerance = 1e-6)
})

# ---- BB1 -----------------------------------------------------------------------

test_that("the BB1 pair copula is supported by the kernels and equals the R functions", {
  set.seed(23)
  for (rot in c(0, 90, 180, 270)) {
    for (par in list(c(0.5, 1.3), c(2, 1), c(1, 2.5), c(0.1, 1.05))) {
      p <- list(family = "bb1", rotation = rot, parameters = par)
      u <- runif(200)
      expect_false(is.null(dv_encode(list(p, p))))
      expect_equal(dv_loglik(u, list(p, p)), dv_loglik_r(u, list(p, p)), tolerance = 1e-6, info = paste(rot, par[1], par[2]))
    }
  }
})

test_that("the gradient with respect to the free BB1 parameter equals finite differences", {
  set.seed(24)
  for (m in bb1_models()) {
    m@K <- 10L
    u <- tssim(m, 600)
    lay <- dv_layout(m)
    z0 <- dv_start(m, NULL, u, list())$z0
    z <- z0 + rnorm(length(z0), 0, 0.03)
    g <- dv_negll_grad(z, lay, 10L, u)
    expect_false(is.null(g))
    fd <- vapply(seq_along(z), function(k) {
      h <- 1e-6
      zp <- z; zp[k] <- zp[k] + h
      zm <- z; zm[k] <- zm[k] - h
      (dv_negll(zp, lay, 10L, u) - dv_negll(zm, lay, 10L, u)) / (2 * h)
    }, 1)
    expect_lt(max(abs(g - fd) / pmax(abs(fd), 1)), 1e-5)
  }
})

test_that("a model with a BB1 lag is fitted with the kernel and the exact gradient and agrees with the R fit", {
  skip_on_cran()
  set.seed(25)
  truth <- sdvinecop(arfactcop(lambda = 0.9, w = 0.5), pairs = list(pc("bb1", aux = 1.2)), tail = pc("joe"), maxlag = 10, tautol = 1e-4)
  u <- tssim(truth, 700)
  spec <- sdvinecop(arfactcop(order = 1), pairs = list(pc("bb1", aux = 1.1)), tail = pc("joe"), maxlag = 10, tautol = 1e-4)
  f1 <- tsfit(u, spec, pseudo = FALSE)
  old <- options(tscop.cpp = FALSE)
  on.exit(options(old))
  f2 <- tsfit(u, spec, pseudo = FALSE)
  expect_equal(unname(coef(f1)), unname(coef(f2)), tolerance = 1e-3)
  expect_equal(as.numeric(logLik(f1)), as.numeric(logLik(f2)), tolerance = 1e-6)
})

# ---- absolute spherical t (ast) -------------------------------------------------

test_that("the ast likelihood kernel equals the R likelihood for plain, IVT and mixed models", {
  set.seed(31)
  for (nm in names(ast_models())) {
    m <- ast_models()[[nm]]
    m@K <- 12L
    pr <- dv_resolve(m)
    u <- tssim(m, 500)
    expect_null(dv_encode(pr))                      # the other kernels do not take ast
    expect_false(is.null(dv_encode(pr, ast = TRUE)))
    expect_equal(dv_loglik(u, pr), dv_loglik_r(u, pr), tolerance = 1e-8, info = nm)
  }
})

test_that("rotated ast pairs fall back to R", {
  set.seed(32)
  m <- sdvinecop(armacop(ar = c(0.5, -0.3)), pairs = list(pc("ast"), pc("ast")), tail = pc("gauss"), maxlag = 6)
  m@K <- 6L
  pr <- dv_resolve(m)
  expect_equal(pr[[2]]$rotation, 90)                # a negative tau at lag 2: a rotated ast pair
  expect_null(dv_encode(pr, ast = TRUE))
  u <- tssim(m, 300)
  expect_equal(dv_loglik(u, pr), dv_loglik_r(u, pr))
  lay <- dv_layout(m)
  z <- dv_start(m, NULL, u, list())$z0
  expect_null(dv_negll_grad(z, lay, 6L, u))         # no exact gradient: central differences are used
  expect_true(all(is.finite(dv_negll_gr(z, lay, 6L, u))))
})

test_that("the exact ast gradient equals large-step differences and differences of the R implementation", {
  set.seed(35)
  cfd <- function(f, z, h) vapply(seq_along(z), function(k) {
    zp <- z; zp[k] <- zp[k] + h
    zm <- z; zm[k] <- zm[k] - h
    (f(zp) - f(zm)) / (2 * h)
  }, 1)
  for (nm in names(ast_models())) {
    m <- ast_models()[[nm]]
    m@K <- 10L
    u <- tssim(m, 500)
    lay <- dv_layout(m)
    z0 <- dv_start(m, NULL, u, list())$z0
    z <- z0 + rnorm(length(z0), 0, 0.03)
    g <- dv_negll_grad(z, lay, 10L, u)
    expect_false(is.null(g), info = nm)
    # Richardson-extrapolated differences with a large step: the compiled value carries the 2e-10
    # error of the level maps, which makes small steps noisy
    f <- function(zz) dv_negll(zz, lay, 10L, u)
    d <- (4 * cfd(f, z, 1e-3) - cfd(f, z, 2e-3)) / 3
    # (models with a lag of nu below 1 are limited to about 1e-5 by R's own qt, accurate to 3e-7 there)
    expect_lt(max(abs(g - d) / pmax(abs(d), 1)), 1e-4)
  }
  # against the R implementation, which uses the exact qt and pt
  m <- ast_models()$ivt
  m@K <- 10L
  u <- tssim(m, 500)
  lay <- dv_layout(m)
  z <- dv_start(m, NULL, u, list())$z0 + rnorm(4, 0, 0.03)
  g <- dv_negll_grad(z, lay, 10L, u)
  old <- options(tscop.cpp = FALSE)
  on.exit(options(old))
  fr <- function(zz) dv_negll(zz, lay, 10L, u)
  dr <- (4 * cfd(fr, z, 1e-4) - cfd(fr, z, 2e-4)) / 3
  expect_lt(max(abs(g - dr) / pmax(abs(dr), 1)), 1e-5)
})

test_that("a model with ast pairs is fitted with the kernel and agrees with the R fit", {
  skip_on_cran()
  set.seed(33)
  kp <- arfactcop(lambda = 0.9, w = 0.5)
  truth <- sdvinecop(kp, tail = pc("ivt", base = "ast"), vt = c(0.45, 0.55), maxlag = 12, tautol = 1e-3)
  u <- tssim(truth, 800)
  spec <- sdvinecop(arfactcop(order = 1), tail = pc("ivt", base = "ast"), maxlag = 12, tautol = 1e-3)
  f1 <- tsfit(u, spec, pseudo = FALSE)
  old <- options(tscop.cpp = FALSE)
  on.exit(options(old))
  f2 <- tsfit(u, spec, pseudo = FALSE)
  expect_equal(unname(coef(f1)), unname(coef(f2)), tolerance = 1e-3)
  expect_equal(as.numeric(logLik(f1)), as.numeric(logLik(f2)), tolerance = 1e-6)
})

# ---- t copula -------------------------------------------------------------------

test_that("the t copula kernel equals the R likelihood for runs of t levels, negative correlation and mixtures", {
  set.seed(41)
  for (nm in names(t_models())) {
    m <- t_models()[[nm]]
    m@K <- 10L
    pr <- dv_resolve(m)
    expect_false(is.null(dv_encode(pr, ast = TRUE)), info = nm)
    expect_equal(dv_loglik(tssim(m, 500) -> u, pr), dv_loglik_r(u, pr), tolerance = 1e-8, info = nm)
  }
})

test_that("the t copula gradient with the free degrees of freedom equals Richardson differences", {
  set.seed(42)
  cfd <- function(f, z, h) vapply(seq_along(z), function(k) {
    zp <- z; zp[k] <- zp[k] + h
    zm <- z; zm[k] <- zm[k] - h
    (f(zp) - f(zm)) / (2 * h)
  }, 1)
  for (nm in names(t_models())) {
    m <- t_models()[[nm]]
    m@K <- 10L
    u <- tssim(m, 500)
    lay <- dv_layout(m)
    z0 <- dv_start(m, NULL, u, list())$z0
    z <- z0 + rnorm(length(z0), 0, 0.03)
    g <- dv_negll_grad(z, lay, 10L, u)
    expect_false(is.null(g), info = nm)
    f <- function(zz) dv_negll(zz, lay, 10L, u)
    d <- (4 * cfd(f, z, 1e-4) - cfd(f, z, 2e-4)) / 3
    expect_lt(max(abs(g - d) / pmax(abs(d), 1)), 1e-5)
  }
})

test_that("a model with t pairs is fitted with the kernel and the exact gradient and agrees with the R fit", {
  skip_on_cran()
  set.seed(43)
  truth <- sdvinecop(arfactcop(lambda = 0.9, w = 0.5), pairs = list(pc("t", aux = 6)), tail = pc("joe"), maxlag = 10, tautol = 1e-3)
  u <- tssim(truth, 700)
  spec <- sdvinecop(arfactcop(order = 1), pairs = list(pc("t", aux = 8)), tail = pc("joe"), maxlag = 10, tautol = 1e-3)
  f1 <- tsfit(u, spec, pseudo = FALSE)
  old <- options(tscop.cpp = FALSE)
  on.exit(options(old))
  f2 <- tsfit(u, spec, pseudo = FALSE)
  expect_equal(unname(coef(f1)), unname(coef(f2)), tolerance = 1e-3)
  expect_equal(as.numeric(logLik(f1)), as.numeric(logLik(f2)), tolerance = 1e-6)
})
