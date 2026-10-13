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
