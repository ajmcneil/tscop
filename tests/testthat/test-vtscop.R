test_that("vtscop validates its base and v-transform", {
  expect_s4_class(vtscop(arfactcop(lambda = 0.9, w = 0.6), udp::v2p(0.45, 0.8)), "vtscop")
  expect_error(vtscop(armacop(ar = 0.5), udp::vlinear(0.5)), "base")
  expect_error(vtscop(arfactcop(order = 1), "x"), "v-transform")
  expect_error(vtscop(arfactcop(order = 1), udp::vlinear(0.5), vtfix = "kappa"), "vtfix")
  expect_error(vtscop(dvinecop(armacop(ar = 0.5), pairs = list(pc("gumbel"))), udp::vlinear(0.5)), "arfactcop")
})

test_that("log-likelihood is the base log-likelihood at V(u)", {
  base <- arfactcop(lambda = 0.9, w = 0.6)
  for (vt in list(udp::vsymmetric(), udp::vlinear(0.4), udp::v2p(0.45, 0.8),
                  udp::v2b(0.45, 0.8), udp::v3p(0.45, 0.8, 1.2), udp::v3b(0.45, 0.8, 1.2))) {
    m <- vtscop(base, vt)
    u <- tssim(m, 300)
    v <- udp::udptrans(vt, u)
    expect_equal(vt_loglik(m, u), vt_base_loglik(base, v))
  }
})

test_that("simulation: V(U) follows the base law and U is uniform", {
  set.seed(11)
  m <- vtscop(arfactcop(lambda = 0.8, w = 0.7), udp::v2p(0.4, 1.5))
  u <- tssim(m, 20000)
  expect_gt(stats::ks.test(u, "punif")$p.value, 0.001)
  v <- udp::udptrans(m@vt, u)
  expect_gt(stats::ks.test(v, "punif")$p.value, 0.001)
  rho1 <- 0.7 * 0.8
  expect_equal(cor(qnorm(v[-1]), qnorm(v[-length(v)])), rho1, tolerance = 0.03)
})

test_that("tsfit recovers parameters and gives standard errors", {
  skip_on_cran()
  set.seed(2)
  m <- vtscop(arfactcop(lambda = 0.9, w = 0.6), udp::v2p(0.45, 0.8))
  u <- tssim(m, 3000)
  f <- tsfit(u, vtscop(arfactcop(order = 1), udp::v2p(0.5, 1)), pseudo = FALSE, se = "hessian")
  expect_s4_class(f, "vtscopfit")
  expect_equal(unname(coef(f)), c(0.9, 0.6, 0.45, 0.8), tolerance = 0.15)
  expect_true(all(is.finite(f@se)) && all(f@se > 0))
  expect_s3_class(logLik(f), "logLik")
  # fixed fulcrum
  g <- tsfit(u, vtscop(arfactcop(order = 1), udp::v2p(0.45, 1), vtfix = "fulcrum"), pseudo = FALSE)
  expect_equal(names(coef(g)), c("lambda1", "w1", "kappa"))
  expect_equal(g@model@vt@pars[["delta"]], 0.45)
})

test_that("tsfit works with a d-vine base and bootstrap standard errors", {
  skip_on_cran()
  set.seed(3)
  m <- vtscop(dvinecop(arfactcop(lambda = 0.9, w = 0.6), pairs = list(pc("gumbel"))), udp::vlinear(0.5))
  u <- tssim(m, 800)
  f <- tsfit(u, vtscop(dvinecop(arfactcop(order = 1), pairs = list(pc("gumbel"))), udp::vlinear(0.5)),
             pseudo = FALSE, se = "bootstrap", B = 10)
  expect_equal(f@se_method, "bootstrap")
  expect_equal(names(coef(f)), c("lambda1", "w1", "fulcrum"))
  expect_true(all(is.finite(f@se)))
})
