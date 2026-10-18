test_that("vtscop validates its base and v-transform", {
  expect_s4_class(vtscop(arfactcop(lambda = 0.9, w = 0.6), udp::v2p(0.45, 0.8)), "vtscop")
  expect_error(vtscop(armacop(ar = 0.5), udp::vlinear(0.5)), "base")
  expect_error(vtscop(arfactcop(order = 1), "x"), "v-transform")
  expect_error(vtscop(arfactcop(order = 1), udp::vlinear(0.5), vtfix = "kappa"), "vtfix")
  expect_error(vtscop(sdvinecop(armacop(ar = 0.5), pairs = list(pc("gumbel"))), udp::vlinear(0.5)), "arfactcop")
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
  g <- tsfit(u, vtscop(arfactcop(order = 1), udp::v2p(0.45, 1), vtfix = "delta"), pseudo = FALSE)
  expect_equal(names(coef(g)), c("lambda1", "w1", "kappa"))
  expect_equal(g@model@vt@pars[["delta"]], 0.45)
})

test_that("tsfit works with a d-vine base and bootstrap standard errors", {
  skip_on_cran()
  set.seed(3)
  m <- vtscop(sdvinecop(arfactcop(lambda = 0.9, w = 0.6), pairs = list(pc("gumbel"))), udp::vlinear(0.5))
  u <- tssim(m, 800)
  f <- tsfit(u, vtscop(sdvinecop(arfactcop(order = 1), pairs = list(pc("gumbel"))), udp::vlinear(0.5)),
             pseudo = FALSE, se = "bootstrap", B = 10)
  expect_equal(f@se_method, "bootstrap")
  expect_equal(names(coef(f)), c("lambda1", "w1", "delta"))
  expect_true(all(is.finite(f@se)))
})

test_that("predict: df, qf and dens are consistent", {
  set.seed(21)
  x <- c(0.05, 0.2, 0.5, 0.6, 0.8, 0.95)
  e <- 1e-5
  for (m in list(vtscop(arfactcop(lambda = 0.9, w = 0.6), udp::vlinear(0.4)),
                 vtscop(arfactcop(lambda = 0.9, w = 0.6), udp::v2p(0.45, 0.8)),
                 vtscop(arfactcop(lambda = 0.9, w = 0.6), udp::v3b(0.45, 0.8, 1.3)),
                 vtscop(sdvinecop(arfactcop(lambda = 0.9, w = 0.6), pairs = list(pc("gumbel"))),
                        udp::v2p(0.45, 0.8)))) {
    u <- tssim(m, 100)
    df <- predict(m, u, type = "df", x = x)
    expect_true(all(diff(df) > 0))
    expect_equal(predict(m, u, type = "qf", x = df), x, tolerance = 1e-6)
    num <- (predict(m, u, type = "df", x = x + e) - predict(m, u, type = "df", x = x - e)) / (2 * e)
    expect_equal(num, predict(m, u, type = "dens", x = x), tolerance = 1e-3)
    expect_equal(predict(m, u, type = "df", x = x, nodes = 150), df, tolerance = 1e-4)
  }
  expect_error(predict(vtscop(arfactcop(lambda = 0.9, w = 0.6), udp::vlinear(0.4)), u, x = 1.2), "inside")
  md <- vtscop(sdvinecop(arfactcop(lambda = 0.9, w = 0.6), pairs = list(pc("gumbel"))), udp::vlinear(0.4))
  expect_error(predict(md, tssim(md, 50), h = 2, x = 0.5), "h = 1")
})

test_that("residuals and predictive PITs are uniform", {
  skip_on_cran()
  set.seed(22)
  m <- vtscop(arfactcop(lambda = 0.9, w = 0.6), udp::v2p(0.45, 0.8))
  U <- vt_simulate(m, 30, 1500)
  pit <- vapply(seq_len(1500), function(j) predict(m, U[1:29, j], type = "df", x = U[30, j]), numeric(1))
  expect_gt(stats::ks.test(pit, "punif")$p.value, 0.001)
  r <- resid(m, U[, 1])
  # predict() tabulates the integral over the base law (accuracy about 1e-5), resid() integrates by quadrature at each point
  expect_equal(r[30], pit[1], tolerance = 1e-5)
  u <- tssim(m, 3000)
  expect_gt(stats::ks.test(resid(m, u), "punif")$p.value, 0.001)
  expect_equal(resid(m, u, scale = "normal"), qnorm(resid(m, u)))
  # d-vine base: linear closed form agrees with the per-step loop result's definition
  md <- vtscop(sdvinecop(arfactcop(lambda = 0.9, w = 0.6), pairs = list(pc("gumbel"))), udp::vlinear(0.4))
  ud <- tssim(md, 1000)
  rd <- resid(md, ud)
  expect_gt(stats::ks.test(rd, "punif")$p.value, 0.001)
  expect_equal(rd[500], predict(md, ud[1:499], type = "df", x = ud[500]), tolerance = 1e-8)
  # fitted object
  f <- tsfit(u, vtscop(arfactcop(order = 1), udp::v2p(0.5, 1)), pseudo = FALSE)
  expect_equal(resid(f), resid(f@model, u))
})

test_that("fast predictive df and quantile functions agree with the quadrature reference", {
  set.seed(7)
  u <- tssim(vtscop(arfactcop(lambda = 0.9, w = 0.6), udp::v2p(0.45, 0.8)), 300)
  xg <- c(0.01, 0.05, 0.2, 0.5, 0.8, 0.95, 0.99)
  ref_df <- function(m, w, x) {
    pd <- vt_predictive_dist(m@base, vt_transform_data(m@vt, w), 1L)
    vt_df(m@vt, x, pd$Fy, pd$Qy, 128)
  }
  for (m in list(vtscop(arfactcop(lambda = 0.9, w = 0.6), udp::v2p(0.45, 0.8)),
                 vtscop(arfactcop(lambda = 0.9, w = 0.6), udp::v3b(0.45, 0.8, 1.3)),
                 vtscop(sdvinecop(arfactcop(lambda = 0.9, w = 0.6), pairs = list(pc("gumbel"))), udp::v2p(0.4, 1.2)))) {
    expect_lt(max(abs(predict(m, u, type = "df", x = xg) - ref_df(m, u, xg))), 5e-5)
    q <- predict(m, u, type = "qf", x = xg)
    expect_true(all(diff(q) > 0))
    expect_lt(max(abs(predict(m, u, type = "df", x = q) - xg)), 1e-9)
  }
  # linear v-transforms: explicit quantile function inverts the df
  ml <- vtscop(sdvinecop(arfactcop(lambda = 0.9, w = 0.6), pairs = list(pc("gumbel"))), udp::vlinear(0.4))
  expect_lt(max(abs(predict(ml, u, type = "df", x = predict(ml, u, type = "qf", x = xg)) - xg)), 1e-9)
})
