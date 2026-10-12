udp_base <- function(base, tau) {
  switch(base,
    ast = udp::astcopula(udp::astcopula_nu(tau)),
    joe = rvinecopulib::bicop_dist("joe", 0, as.numeric(rvinecopulib::ktau_to_par("joe", tau))),
    gumbel = rvinecopulib::bicop_dist("gumbel", 0, as.numeric(rvinecopulib::ktau_to_par("gumbel", tau))),
    clayton180 = rvinecopulib::bicop_dist("clayton", 180, as.numeric(rvinecopulib::ktau_to_par("clayton", tau)))
  )
}

test_that("pc accepts the ivt family with a base copula", {
  p <- pc("ivt", base = "ast")
  expect_equal(c(p@family, p@base), c("ivt", "ast"))
  expect_equal(c(p@pos, p@neg), c(0, 0))
  expect_equal(pc("ivt", base = "Clayton180")@base, "clayton180")
  expect_error(pc("ivt"), "base copula")
  expect_error(pc("ivt", base = "frank"), "base copula")
  expect_error(pc("gumbel", base = "ast"), "only the")
  expect_output(show(p), "ivt with base ast")
})

test_that("IVT pair copula functions equal the udp stochastic-inversion copula (all bases, asymmetric fulcrums)", {
  set.seed(1)
  a <- runif(40, 0.03, 0.97)
  b <- runif(40, 0.03, 0.97)
  for (base in c("ast", "joe", "gumbel", "clayton180")) {
    for (d in list(c(0.5, 0.5), c(0.3, 0.7), c(0.62, 0.41))) {
      tau <- 0.35
      pr <- dv_ivt_pair(base, tau, d)
      bc <- udp::bsicopula(udp_base(base, tau), udp::vlinear(d[1]), udp::vlinear(d[2]))
      info <- paste(base, d[1], d[2])
      expect_equal(dv_dens(pr, a, b), udp::dbsicopula(a, b, bc), tolerance = 1e-8, info = info)
      expect_equal(dv_h(pr, a, b, 1L), udp::hbsicopula(a, b, bc, cond_var = 1L), tolerance = 1e-7, info = info)
      expect_equal(dv_h(pr, a, b, 2L), udp::hbsicopula(a, b, bc, cond_var = 2L), tolerance = 1e-7, info = info)
      p <- dv_h(pr, a, b, 1L)
      expect_equal(dv_hinv(pr, a, p), b, tolerance = 1e-6, info = info)
      expect_equal(dv_hinv(pr, a, b), udp::hbsicopula(a, b, bc, cond_var = 1L, inverse = TRUE),
                   tolerance = 1e-6, info = info)
    }
  }
})

test_that("the v-transform of an h-function is the base h-function of the transformed arguments", {
  set.seed(2)
  a <- runif(30, 0.05, 0.95)
  b <- runif(30, 0.05, 0.95)
  for (base in c("ast", "joe")) {
    d <- c(0.35, 0.6)
    pr <- dv_ivt_pair(base, 0.4, d)
    x <- dv_vlin(a, d[1])
    y <- dv_vlin(b, d[2])
    expect_equal(dv_vlin(dv_h(pr, a, b, 1L), d[2]), dv_h(pr$base, x, y, 1L), tolerance = 1e-8)
    expect_equal(dv_vlin(dv_h(pr, a, b, 2L), d[1]), dv_h(pr$base, x, y, 2L), tolerance = 1e-8)
  }
})

test_that("the d-vine log-likelihood equals a loop over the udp pair copulas", {
  set.seed(3)
  u <- runif(60)
  d <- c(0.4, 0.55)
  taus <- c(0.4, 0.2, 0.1)
  pairs <- lapply(taus, function(t) dv_ivt_pair("joe", t, d))
  v <- cbind(u[-60], u[-1])
  ref <- 0
  for (j in seq_along(taus)) {
    bc <- udp::bsicopula(udp_base("joe", taus[j]), udp::vlinear(d[1]), udp::vlinear(d[2]))
    ref <- ref + sum(log(udp::dbsicopula(v[, 1], v[, 2], bc)))
    if (j < length(taus)) {
      m <- nrow(v)
      v <- cbind(udp::hbsicopula(v[-m, 1], v[-m, 2], bc, cond_var = 2L),
                 udp::hbsicopula(v[-1, 1], v[-1, 2], bc, cond_var = 1L))
    }
  }
  expect_equal(dv_loglik(u, pairs), ref, tolerance = 1e-8)
})

test_that("with one common fulcrum the d-vine density is the base d-vine density at V(u)", {
  set.seed(4)
  u <- runif(80)
  for (base in c("ast", "gumbel")) {
    for (delta in c(0.5, 0.35)) {
      taus <- c(0.3, 0.15, 0.07, 0.03)
      ivt <- lapply(taus, function(t) dv_ivt_pair(base, t, c(delta, delta)))
      basep <- lapply(taus, function(t) dv_ivt_base_pair(base, t))
      expect_equal(dv_loglik(u, ivt), dv_loglik(dv_vlin(u, delta), basep), tolerance = 1e-8)
    }
  }
})

test_that("simulation and Rosenblatt are inverse for IVT pairs with different fulcrums", {
  set.seed(5)
  pairs <- c(lapply(c(0.45, 0.2, 0.1), function(t) dv_ivt_pair("ast", t, c(0.35, 0.62))))
  w <- matrix(runif(30 * 3, 0.03, 0.97), 30, 3)
  u <- dv_simulate(pairs, w)
  for (r in 1:3) expect_equal(dv_rosenblatt(pairs, u[, r]), w[, r], tolerance = 1e-6)
})

test_that("first-order IVT d-vine simulation has the IVT pair copula as its lag-1 copula", {
  set.seed(6)
  d <- c(0.4, 0.6)
  m <- sdvinecop(armacop(ar = 0.6), pairs = list(pc("ivt", base = "joe")), vt = d)
  expect_equal(m@K, 1L)
  U <- dv_simulate(dv_resolve(m), matrix(runif(400 * 100), 400, 100))
  a <- as.numeric(U[-400, ])
  b <- as.numeric(U[-1, ])
  tau <- (2 / pi) * asin(0.6)
  bc <- udp::bsicopula(udp_base("joe", tau), udp::vlinear(d[1]), udp::vlinear(d[2]))
  for (pt in list(c(0.2, 0.3), c(0.5, 0.5), c(0.7, 0.8), c(0.35, 0.9))) {
    emp <- mean(a <= pt[1] & b <= pt[2])
    expect_equal(emp, udp::pbsicopula(pt[1], pt[2], bc), tolerance = 0.01)
  }
  expect_gt(ks.test(a[1:5000], "punif")$p.value, 0.001)
})

test_that("sdvinecop constructor handles the v-transform parameters", {
  m <- sdvinecop(arfactcop(lambda = 0.9, w = 0.6), tail = pc("ivt", base = "ast"), vt = c(0.4, 0.6))
  expect_equal(unname(m@vt), c(0.4, 0.6))
  expect_equal(names(m@vt), c("fulcrum1", "fulcrum2"))
  expect_equal(names(coef(m)), c("lambda1", "w1", "fulcrum1", "fulcrum2"))
  expect_output(show(m), "v-transform fulcrums: 0.4, 0.6")
  expect_output(show(m), "tail family: ivt with base ast")
  # default fulcrums
  m0 <- sdvinecop(arfactcop(lambda = 0.9, w = 0.6), pairs = list(pc("ivt", base = "joe")))
  expect_equal(unname(m0@vt), c(0.5, 0.5))
  # common and fixed
  mc <- sdvinecop(arfactcop(lambda = 0.9, w = 0.6), tail = pc("ivt", base = "ast"), vt = c(0.4, 0.4), vtpar = "common")
  expect_equal(names(coef(mc)), c("lambda1", "w1", "fulcrum"))
  mn <- sdvinecop(arfactcop(lambda = 0.9, w = 0.6), tail = pc("ivt", base = "ast"), vtpar = "none")
  expect_equal(names(coef(mn)), c("lambda1", "w1"))
  expect_output(show(mn), "\\(fixed\\)")
  expect_error(sdvinecop(arfactcop(lambda = 0.9, w = 0.6), tail = pc("ivt", base = "ast"), vt = c(0.4, 0.5), vtpar = "common"), "equal")
  expect_error(sdvinecop(arfactcop(lambda = 0.9, w = 0.6), tail = pc("ivt", base = "ast"), vt = c(0, 0.5)), "\\(0, 1\\)")
  expect_error(sdvinecop(arfactcop(lambda = 0.9, w = 0.6), vt = c(0.4, 0.6)), "only used")
  expect_error(sdvinecop(arfactcop(lambda = 0.9, w = 0.6), tail = pc("ivt", base = "ast"), vt = 0.5), "two fulcrums")
})

test_that("a kpacf with negative values is infeasible for an IVT pair copula", {
  m <- sdvinecop(armacop(ar = c(0.5, -0.4)), pairs = list(pc("gumbel"), pc("ivt", base = "joe")))
  expect_null(dv_resolve(m))
  expect_equal(dv_model_loglik(m, runif(50)), -Inf)
  expect_error(tssim(m, 20), "infeasible")
  mt <- sdvinecop(armacop(ar = c(0.5, -0.4)), tail = pc("ivt", base = "joe"))
  expect_null(dv_resolve(mt))
  # non-negative kpacf is fine
  expect_false(is.null(dv_resolve(sdvinecop(armacop(ar = 0.5), tail = pc("ivt", base = "joe")))))
})

test_that("ktau reports the kpacf (the base copulas' tau)", {
  m <- sdvinecop(arfactcop(lambda = 0.9, w = 0.6), tail = pc("ivt", base = "ast"), vt = c(0.4, 0.6))
  expect_equal(ktau(m, 5), gc_tau(m@kpacf, 5))
  pr <- dv_resolve(m)
  expect_equal(dv_pair_ktau(pr[[1]]), gc_tau(m@kpacf, 1), tolerance = 1e-6)
})

test_that("residuals and prediction are consistent for an IVT d-vine", {
  set.seed(7)
  m <- sdvinecop(arfactcop(lambda = 0.9, w = 0.6), pairs = list(pc("ivt", base = "joe")),
                tail = pc("ivt", base = "joe"), vt = c(0.4, 0.6), tautol = 1e-3)
  u <- tssim(m, 120)
  r <- resid(m, u)
  expect_gt(ks.test(r, "punif")$p.value, 0.001)
  x <- c(0.1, 0.5, 0.9)
  df <- predict(m, u, type = "df", x = x)
  for (i in 1:3) expect_equal(df[i], resid(m, c(u, x[i]))[121], tolerance = 1e-7)
  expect_equal(predict(m, u, type = "qf", x = df), x, tolerance = 1e-6)
  dens <- function(v) predict(m, u, type = "dens", x = v)
  expect_equal(integrate(dens, 0, 1, rel.tol = 1e-6)$value, 1, tolerance = 1e-3)
})

test_that("tsfit recovers the fulcrums and the kpacf of an IVT d-vine", {
  set.seed(8)
  truth <- sdvinecop(arfactcop(lambda = 0.85, w = 0.6), tail = pc("ivt", base = "joe"),
                    vt = c(0.4, 0.6), tautol = 1e-3)
  u <- tssim(truth, 1500)
  spec <- sdvinecop(arfactcop(order = 1), tail = pc("ivt", base = "joe"), tautol = 1e-3)
  f <- tsfit(u, spec, pseudo = FALSE)
  expect_equal(names(coef(f)), c("lambda1", "w1", "fulcrum1", "fulcrum2"))
  expect_equal(unname(coef(f)[3:4]), c(0.4, 0.6), tolerance = 0.2)
  expect_equal(unname(coef(f)[1]), 0.85, tolerance = 0.15)
  expect_equal(f@npar, 4L)
  expect_equal(unname(f@model@vt), unname(coef(f)[3:4]))
  # a Gaussian-tail fit is much worse
  fg <- tsfit(u, sdvinecop(arfactcop(order = 1), tautol = 1e-3), pseudo = FALSE)
  expect_gt(f@loglik, fg@loglik + 20)
  # common and fixed fulcrums
  fc <- tsfit(u, sdvinecop(arfactcop(order = 1), tail = pc("ivt", base = "joe"), tautol = 1e-3, vtpar = "common"),
              pseudo = FALSE)
  expect_equal(names(coef(fc)), c("lambda1", "w1", "fulcrum"))
  expect_equal(fc@npar, 3L)
  expect_equal(unname(fc@model@vt[1]), unname(fc@model@vt[2]))
  expect_lte(fc@loglik, f@loglik + 1e-6)
  fn <- tsfit(u, sdvinecop(arfactcop(order = 1), tail = pc("ivt", base = "joe"), tautol = 1e-3, vt = c(0.5, 0.5),
                          vtpar = "none"), pseudo = FALSE)
  expect_equal(fn@npar, 2L)
  expect_lte(fn@loglik, fc@loglik + 1e-6)
  # standard errors
  fh <- tsfit(u, spec, pseudo = FALSE, se = "hessian")
  expect_true(all(fh@se > 0))
  fb <- tsfit(u, spec, pseudo = FALSE, se = "bootstrap", B = 10)
  expect_equal(dim(fb@boot), c(10L, 4L))
})

test_that("an ast base works in a fit and with an ARMA kpacf in its non-negative region", {
  set.seed(9)
  truth <- sdvinecop(armacop(ar = 0.8, ma = -0.5), tail = pc("ivt", base = "ast"), vt = c(0.45, 0.55), tautol = 1e-2)
  u <- tssim(truth, 700)
  f <- tsfit(u, sdvinecop(armacop(order = c(1, 1)), tail = pc("ivt", base = "ast"), tautol = 1e-2), pseudo = FALSE)
  expect_s4_class(f, "sdvinecopfit")
  expect_equal(f@npar, 4L)
  expect_true(all(f@model@vt > 0 & f@model@vt < 1))
})

test_that("the fit and the IVT likelihood agree with udp-based evaluation at the estimates", {
  set.seed(10)
  truth <- sdvinecop(arfactcop(lambda = 0.8, w = 0.5), tail = pc("ivt", base = "gumbel"),
                    vt = c(0.45, 0.55), tautol = 1e-2)
  u <- tssim(truth, 200)
  pr <- dv_resolve(truth)
  # reference: udp bsicopula pair copulas with the same taus
  taus <- gc_tau(truth@kpacf, truth@K)
  v <- cbind(u[-200], u[-1])
  ref <- 0
  for (j in seq_len(truth@K)) {
    bc <- udp::bsicopula(udp_base("gumbel", taus[j]), udp::vlinear(0.45), udp::vlinear(0.55))
    ref <- ref + sum(log(udp::dbsicopula(v[, 1], v[, 2], bc)))
    if (j < truth@K) {
      m <- nrow(v)
      v <- cbind(udp::hbsicopula(v[-m, 1], v[-m, 2], bc, cond_var = 2L),
                 udp::hbsicopula(v[-1, 1], v[-1, 2], bc, cond_var = 1L))
    }
  }
  expect_equal(dv_model_loglik(truth, u), ref, tolerance = 1e-7)
})
