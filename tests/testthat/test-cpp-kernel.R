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

test_that("unsupported families use the R implementation and the option switches the kernel off", {
  kp <- arfactcop(lambda = 0.9, w = 0.5)
  m <- sdvinecop(kp, pairs = list(pc("t", aux = 6), pc("gumbel")), tail = "gauss", maxlag = 6)
  m@K <- 6L
  pr <- dv_resolve(m)
  expect_null(dv_encode(pr))
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
  m <- sdvinecop(arfactcop(lambda = 0.9, w = 0.5), pairs = list(pc("t", aux = 6)), tail = "gauss", maxlag = 6)
  m@K <- 6L
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

test_that("vt_negll_grad declines a Gaussian base and vt_negll_gr then uses differences", {
  set.seed(10)
  m <- vtscop(arfactcop(lambda = 0.9, w = 0.5), udp::v2p(0.45, 0.8))
  u <- tssim(m, 300)
  L <- vt_layout(m)
  z <- vt_start(m, NULL, u, list())$z0
  expect_null(vt_negll_grad(z, L, NA_integer_, u))
  expect_true(all(is.finite(vt_negll_gr(z, L, NA_integer_, u))))
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
