aicc_formula <- function(f) {
  ll <- logLik(f)
  k <- attr(ll, "df")
  n <- attr(ll, "nobs")
  -2 * as.numeric(ll) + 2 * k + 2 * k * (k + 1) / (n - k - 1)
}

test_that("AICc is the corrected AIC for every fitted model class", {
  set.seed(1)
  kp <- arfactcop(lambda = 0.9, w = 0.5)
  u <- tssim(sdvinecop(kp, pairs = list(pc("gumbel")), maxlag = 10, tautol = 1e-3), 300)
  fits <- list(
    armacop = tsfit(u, armacop(order = c(1, 1)), pseudo = FALSE),
    arfactcop = tsfit(u, arfactcop(order = 1), pseudo = FALSE),
    sdvinecop = tsfit(u, sdvinecop(arfactcop(order = 1), pairs = list(pc("gumbel")), maxlag = 10, tautol = 1e-3), pseudo = FALSE),
    vtscop = tsfit(u, vtscop(arfactcop(order = 1), udp::v2p(0.5, 1)), pseudo = FALSE),
    margfit = tsfit(rnorm(200), marg("st")),
    tsmod = tsfit(tssim(tsmod(marg("st", c(mu = 0, sigma = 1, df = 6)), kp), 300), tsmod("st", arfactcop(order = 1))),
    tsmod_emarg = tsfit(rnorm(300), tsmod(emarg(rnorm(300)), arfactcop(order = 1)))
  )
  for (nm in names(fits)) {
    f <- fits[[nm]]
    expect_equal(AICc(f), aicc_formula(f), info = nm)
    expect_gt(AICc(f), AIC(f))
    expect_equal(AICc(f) - AIC(f), 2 * attr(logLik(f), "df") * (attr(logLik(f), "df") + 1) /
                   (attr(logLik(f), "nobs") - attr(logLik(f), "df") - 1), info = nm)
  }
})

test_that("AICc of several models is a data frame with the call's names, as AIC", {
  set.seed(2)
  u <- tssim(armacop(ar = 0.6, ma = 0.3), 300)
  f1 <- tsfit(u, armacop(order = c(1, 0)), pseudo = FALSE)
  f2 <- tsfit(u, armacop(order = c(1, 1)), pseudo = FALSE)
  r <- AICc(f1, f2)
  expect_s3_class(r, "data.frame")
  expect_equal(names(r), c("df", "AICc"))
  expect_equal(rownames(r), c("f1", "f2"))
  expect_equal(r$AICc, c(AICc(f1), AICc(f2)))
  expect_equal(r$df, c(1, 2))
  expect_equal(rownames(AIC(f1, f2)), rownames(r))
})

test_that("AICc is infinite when the correction is undefined and works for other logLik classes", {
  expect_equal(ic_aicc(-10, 5, 6), Inf)
  expect_equal(ic_aicc(-10, 5, 5), Inf)
  expect_equal(ic_aicc(-10, 2, 20), 20 + 4 + 2 * 2 * 3 / 17)
  fit <- lm(dist ~ speed, cars)
  k <- attr(logLik(fit), "df")
  expect_equal(AICc(fit), AIC(fit) + 2 * k * (k + 1) / (nrow(cars) - k - 1))
})

test_that("automarg ranks the margin families by AICc on request", {
  set.seed(3)
  x <- tssim(marg("st", c(mu = 0, sigma = 1, df = 4)), 40)          # a small sample: AICc differs from AIC
  a <- automarg(x, families = c("gauss", "st", "sst", "gh"), criterion = "AICc")
  expect_equal(a@criterion, "AICc")
  expect_true("AICc" %in% names(a@table))
  ok <- !is.na(a@table$loglik)
  expect_equal(a@table$AICc[ok], ic_aicc(a@table$loglik[ok], a@table$npar[ok], length(x)))
  expect_equal(a@table$family[1], a@table$family[which.min(a@table$AICc)])
  expect_equal(a@table$delta[1], 0)
  expect_error(automarg(x, criterion = "foo"))
})
