# Numerical comparison with the legacy tscopula package (vtscopula with an
# ARMA(1,1) base). An arfactcop with one factor and noise is an ARMA(1,1)
# copula with ar = lambda, ma in [-ar, 0] (see arfact2arma), so the models
# coincide when the ARMA parameters lie in that region.

legacy_to_vtscop <- function(g) {
  cg <- tscopula::coef(g)
  ar <- cg[["ar1"]]
  ma <- cg[["ma1"]]
  rho1 <- (1 + ar * ma) * (ar + ma) / (1 + 2 * ar * ma + ma^2)
  vtscop(arfactcop(lambda = ar, w = rho1 / ar), udp::v2p(cg[["delta"]], cg[["kappa"]]))
}

test_that("vtscop fits and log-likelihoods agree with tscopula::vtscopula", {
  skip_if_no_oracle()
  skip_on_cran()
  set.seed(7)
  u <- tssim(vtscop(arfactcop(lambda = 0.9, w = 0.6), udp::v2p(0.5, 0.8)), 1500)
  g <- tscopula::fit(tscopula::vtscopula(tscopula::armacopula(list(ar = 0.5, ma = -0.2)),
                                         tscopula::V2p(0.5, 1)), u)
  f <- tsfit(u, vtscop(arfactcop(order = 1), udp::v2p(0.5, 1), vtfix = "fulcrum"), pseudo = FALSE)
  # same maximised likelihood and parameters (via the ARMA(1,1) representation)
  expect_equal(as.numeric(logLik(f)), as.numeric(logLik(g)), tolerance = 1e-4)
  a <- coef(arfact2arma(f@model@base))
  cg <- tscopula::coef(g)
  expect_equal(unname(a), unname(cg[c("ar1", "ma1")]), tolerance = 1e-2)
  expect_equal(f@model@vt@pars[["delta"]], cg[["delta"]])
  expect_equal(f@model@vt@pars[["kappa"]], cg[["kappa"]], tolerance = 1e-2)
  # log-likelihood of our model at the legacy estimates equals the legacy value
  mine <- legacy_to_vtscop(g)
  expect_equal(vt_loglik(mine, u), as.numeric(logLik(g)), tolerance = 1e-6)
  # kpacf
  expect_equal(ktau(mine, 6), as.numeric(tscopula::kendall(g@tscopula, 6)), tolerance = 1e-6)
})

test_that("tscopula residuals are the base-process residuals of V(u)", {
  skip_if_no_oracle()
  skip_on_cran()
  # tscopula::resid on a fitted vtscopula returns the normal-scale residuals of
  # the base process at V(u); they differ from the first values only because of
  # its start-up approximation. (tscop's resid is the conditional PIT of u.)
  set.seed(8)
  u <- tssim(vtscop(arfactcop(lambda = 0.9, w = 0.6), udp::v2p(0.5, 0.8)), 1000)
  g <- tscopula::fit(tscopula::vtscopula(tscopula::armacopula(list(ar = 0.5, ma = -0.2)),
                                         tscopula::V2p(0.5, 1)), u)
  mine <- legacy_to_vtscop(g)
  v <- udp::udptrans(mine@vt, u)
  rb <- resid(mine@base, v, scale = "normal")
  old <- as.numeric(tscopula::resid(g))
  expect_equal(rb[-(1:20)], old[-(1:20)], tolerance = 1e-6)
  # and our residuals are uniform for the same model
  expect_gt(stats::ks.test(resid(mine, u), "punif")$p.value, 0.001)
})
