test_that("tssim simulates the factor process exactly from the stationary distribution", {
  set.seed(1)
  g <- arfactcop(lambda = c(0.95, 0.6), w = c(0.3, 0.4))
  X <- qnorm(t(vapply(1:15000, function(i) tssim(g, 3), numeric(3))))
  rho <- afc_acf(g@lambda, g@w, 2)
  expect_equal(unname(cov(X)), toeplitz(c(1, rho)), tolerance = 0.06, ignore_attr = TRUE)
  x <- qnorm(tssim(g, 2e5))
  expect_equal(var(x), 1, tolerance = 0.03)
  expect_equal(as.numeric(acf(x, lag.max = 4, plot = FALSE)$acf[-1]), afc_acf(g@lambda, g@w, 4), tolerance = 0.03)
  set.seed(2)
  a <- tssim(g, 10)
  set.seed(2)
  expect_equal(a, tssim(g, 10))
  expect_length(tssim(g), 1000)
  nn <- arfactcop(lambda = 0.8, w = 1, noise = FALSE)
  expect_equal(var(qnorm(tssim(nn, 1e5))), 1, tolerance = 0.03)
  expect_error(tssim(g, 0), "positive")
  expect_error(tssim(arfactcop(order = 2), 10), "unspecified")
})

test_that("ktau of an arfactcop is the Kendall's tau of its PACF", {
  g <- arfactcop(lambda = c(0.95, 0.6), w = c(0.3, 0.4))
  expect_equal(ktau(g, 6), gc_tau(g, 6))
  expect_length(ktau(g), 20)
  expect_error(ktau(g, 0), "lagmax")
})

test_that("resid gives the Cholesky-standardised normal scores", {
  set.seed(3)
  g <- arfactcop(lambda = c(0.95, 0.6), w = c(0.3, 0.4))
  u <- tssim(g, 50)
  L <- t(chol(toeplitz(c(1, afc_acf(g@lambda, g@w, 49)))))
  ref <- as.numeric(forwardsolve(L, qnorm(u)))
  expect_equal(resid(g, u, scale = "normal"), ref, tolerance = 1e-8)
  expect_equal(resid(g, u), pnorm(ref), tolerance = 1e-8)
  expect_equal(residuals(g, u), resid(g, u))
  u2 <- u
  u2[20] <- NA
  r <- resid(g, u2, scale = "normal")
  expect_true(is.na(r[20]))
  expect_equal(r[1:19], ref[1:19], tolerance = 1e-8)
  ut <- ts(u, start = c(2000, 1), frequency = 12)
  expect_s3_class(resid(g, ut), "ts")
  # iid for data from the model
  uu <- tssim(g, 20000)
  z <- resid(g, uu, scale = "normal")
  expect_lt(abs(mean(z)), 0.03)
  expect_equal(var(z), 1, tolerance = 0.05)
  expect_lt(max(abs(acf(z, lag.max = 8, plot = FALSE)$acf[-1])), 0.03)
})

test_that("predict matches exact Gaussian conditioning and is internally consistent", {
  set.seed(4)
  g <- arfactcop(lambda = c(0.95, 0.6), w = c(0.3, 0.4))
  u <- tssim(g, 30)
  acf <- c(1, afc_acf(g@lambda, g@w, 60))
  pr <- predict(g, u, h = c(1, 2, 5))
  for (i in 1:3) {
    ref <- cond_moments(qnorm(u), acf, pr$h[i])
    expect_equal(c(pr$mean[i], pr$sd[i]), unname(ref), tolerance = 1e-8)
  }
  um <- u
  um[c(10, 29, 30)] <- NA
  prm <- predict(g, um, h = 1)
  expect_equal(c(prm$mean, prm$sd), unname(cond_moments(qnorm(um), acf, 1)), tolerance = 1e-8)
  x <- c(0.05, 0.3, 0.6, 0.95)
  df <- predict(g, u, type = "df", x = x)
  expect_equal(predict(g, u, type = "qf", x = df), x, tolerance = 1e-8)
  dens <- function(v) predict(g, u, type = "dens", x = v)
  expect_equal(integrate(dens, 0, 1, rel.tol = 1e-8)$value, 1, tolerance = 1e-4)
  expect_error(predict(g, u, h = 0), "'h'")
  expect_error(predict(g, u, type = "df"), "'x'")
})

test_that("tsfit recovers the factor model and agrees with armacop for one factor", {
  x <- sim_arfact_scores(2500, c(0.97, 0.6), c(0.3, 0.3), seed = 5)
  u <- pnorm(x)
  f <- tsfit(u, arfactcop(order = 2), pseudo = FALSE)
  expect_s4_class(f, "arfactcopfit")
  expect_equal(names(coef(f)), c("lambda1", "lambda2", "w1", "w2"))
  expect_equal(f@loglik, afc_copula_loglik(u, f@model@lambda, f@model@w, 1 - sum(f@model@w)), tolerance = 1e-8)
  expect_equal(as.numeric(logLik(f)), f@loglik)
  expect_equal(f@npar, 4L)
  expect_equal(AIC(f), -2 * f@loglik + 8)
  # one factor: arfactcop is a sub-model of ARMA(1,1)
  x1 <- sim_arfact_scores(2000, 0.9, 0.6, seed = 6)
  u1 <- pnorm(x1)
  f1 <- tsfit(u1, arfactcop(order = 1), pseudo = FALSE)
  fa <- tsfit(u1, armacop(order = c(1, 1)), pseudo = FALSE)
  expect_lte(f1@loglik, fa@loglik + 1e-6)
  expect_equal(f1@loglik, fa@loglik, tolerance = 2e-3)
  # known model as the start; init as a numeric vector; NA values
  expect_equal(tsfit(u1, f1@model, pseudo = FALSE)@loglik, f1@loglik, tolerance = 1e-6)
  expect_equal(tsfit(u1, arfactcop(order = 1), init = c(0.5, 0.2), pseudo = FALSE)@loglik, f1@loglik, tolerance = 1e-5)
  u1m <- u1
  u1m[c(5, 6, 700)] <- NA
  expect_equal(tsfit(u1m, arfactcop(order = 1), pseudo = FALSE)@nobs, 1997L)
  expect_error(tsfit(u1, arfactcop(order = 1), se = "foo"), "se")
  expect_error(tsfit(u1, arfactcop(order = 1), se = "bootstrap", B = 3), "B")
  expect_output(show(f1), "arfactcopfit")
})

test_that("standard errors for the factor model", {
  x1 <- sim_arfact_scores(1500, 0.9, 0.6, seed = 7)
  u1 <- pnorm(x1)
  h <- tsfit(u1, arfactcop(order = 1), pseudo = FALSE, se = "hessian")
  expect_equal(h@se_method, "hessian")
  expect_true(all(h@se > 0))
  expect_equal(vcov(h), t(vcov(h)))
  b <- tsfit(u1, arfactcop(order = 1), pseudo = FALSE, se = "bootstrap", B = 30)
  expect_equal(b@se_method, "bootstrap")
  expect_equal(dim(b@boot), c(30L, 2L))
  expect_equal(unname(b@se), unname(h@se), tolerance = 0.6)
  r <- tsfit(u1, arfactcop(order = 1), se = TRUE, B = 10)
  expect_equal(dim(r@boot), c(10L, 2L))
})

test_that("methods on the fitted object delegate to the model", {
  x1 <- sim_arfact_scores(400, 0.9, 0.6, seed = 8)
  u1 <- pnorm(x1)
  f <- tsfit(u1, arfactcop(order = 1), pseudo = FALSE)
  expect_equal(ktau(f, 3), ktau(f@model, 3))
  expect_equal(resid(f), resid(f@model, u1))
  expect_equal(residuals(f, scale = "normal"), resid(f@model, u1, scale = "normal"))
  expect_equal(predict(f, h = 1:2), predict(f@model, u1, h = 1:2))
  expect_length(tssim(f), 400)
  expect_length(tssim(f, 20), 20)
  ft <- tsfit(ts(u1, start = c(2000, 1), frequency = 12), arfactcop(order = 1), pseudo = FALSE)
  expect_s3_class(resid(ft), "ts")
})

test_that("residuals agree with the legacy ARMA(1,1) copula once its filter has settled", {
  skip_if_no_oracle()
  skip_on_cran()
  # one factor plus noise is an ARMA(1,1) copula with ar = lambda, -ar <= ma <= 0;
  # tscopula's residuals are only accurate after a start-up stretch
  set.seed(12)
  u <- tssim(arfactcop(lambda = 0.9, w = 0.6), 300)
  g <- tscopula::fit(tscopula::armacopula(list(ar = 0.5, ma = -0.2)), u)
  ar <- tscopula::coef(g)[["ar1"]]
  ma <- tscopula::coef(g)[["ma1"]]
  rho1 <- (1 + ar * ma) * (ar + ma) / (1 + 2 * ar * ma + ma^2)
  mine <- arfactcop(lambda = ar, w = rho1 / ar)
  expect_equal(resid(mine, u, scale = "normal")[21:300], as.numeric(tscopula::resid(g))[21:300],
               tolerance = 1e-6)
})
