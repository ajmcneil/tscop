test_that("free <-> natural maps round-trip across blocks", {
  pars <- list(ar = c(0.5, -0.2), ma = 0.4, sar = 0.3, sma = c(-0.2, 0.1))
  ord <- c(2L, 1L, 1L, 2L)
  z <- pars_to_free(pars)
  e <- free_to_expanded(z, ord, 12)
  ref <- arma_expand(pars$ar, pars$ma, pars$sar, pars$sma, 12)
  expect_equal(e$ar, ref$ar)
  expect_equal(e$ma, ref$ma)
})

test_that("Hannan-Rissanen starts are admissible and near the truth", {
  x <- sim_arma_scores(3000, ar = 0.6, ma = 0.3, seed = 1)
  s <- arma_hr_start(x, c(1L, 1L, 0L, 0L), NA_real_)
  expect_equal(s$ar, 0.6, tolerance = 0.15)
  expect_equal(s$ma, 0.3, tolerance = 0.3)
  expect_true(arma_is_stationary(s$ar))
  expect_true(arma_is_invertible(s$ma))
  expect_length(s$sar, 0)
  # seasonal
  xs <- sim_arma_scores(2000, ar = arma_expand(sar = 0.5, period = 4)$ar, seed = 2)
  ss <- arma_hr_start(xs, c(0L, 0L, 1L, 0L), 4)
  expect_equal(ss$sar, 0.5, tolerance = 0.1)
})

test_that("shrinking makes inadmissible starts admissible", {
  a <- arma_shrink(c(1.2, 0.3))
  expect_true(arma_is_stationary(a))
  expect_equal(arma_shrink(0.5), 0.5)
})

test_that("tsfit recovers AR(1), MA(1) and ARMA(1,1) parameters", {
  x <- sim_arma_scores(3000, ar = 0.6, seed = 3)
  f <- tsfit(pnorm(x), armacop(order = c(1, 0)), pseudo = FALSE)
  expect_s4_class(f, "armacopfit")
  expect_equal(unname(coef(f)), 0.6, tolerance = 0.06)
  expect_equal(f@convergence, 0L)
  x <- sim_arma_scores(3000, ma = 0.5, seed = 4)
  f <- tsfit(pnorm(x), armacop(order = c(0, 1)), pseudo = FALSE)
  expect_equal(unname(coef(f)), 0.5, tolerance = 0.08)
  x <- sim_arma_scores(4000, ar = 0.5, ma = 0.4, seed = 5)
  f <- tsfit(pnorm(x), armacop(order = c(1, 1)), pseudo = FALSE)
  expect_equal(unname(coef(f)), c(0.5, 0.4), tolerance = 0.12)
})

test_that("tsfit recovers a seasonal model", {
  e <- arma_expand(ar = 0.4, sar = 0.5, period = 4)
  x <- sim_arma_scores(3000, ar = e$ar, seed = 6)
  f <- tsfit(pnorm(x), armacop(order = c(1, 0), seasonal = c(1, 0), period = 4),
             pseudo = FALSE)
  expect_equal(unname(coef(f)), c(0.4, 0.5), tolerance = 0.1)
  expect_equal(f@model@period, 4)
})

test_that("the fitted log-likelihood is the copula log-likelihood at the estimate", {
  x <- sim_arma_scores(500, ar = 0.5, ma = 0.3, seed = 7)
  u <- pnorm(x)
  f <- tsfit(u, armacop(order = c(1, 1)), pseudo = FALSE)
  e <- arma_coef_expanded(f@model)
  expect_equal(f@loglik, arma_copula_loglik(u, e$ar, e$ma))
  expect_equal(as.numeric(logLik(f)), f@loglik)
  expect_equal(attr(logLik(f), "df"), 2)
  expect_equal(attr(logLik(f), "nobs"), 500)
  expect_equal(AIC(f), -2 * f@loglik + 4)
})

test_that("the fit is at least as good as the starting values", {
  x <- sim_arma_scores(800, ar = 0.5, ma = 0.4, seed = 8)
  u <- pnorm(x)
  f <- tsfit(u, armacop(order = c(1, 1)), pseudo = FALSE)
  e <- arma_coef_expanded(f@start)
  expect_gte(f@loglik, arma_copula_loglik(u, e$ar, e$ma) - 1e-8)
})

test_that("different starting values reach the same optimum", {
  x <- sim_arma_scores(1000, ar = c(0.5, 0.3), ma = 0.4, seed = 9)
  u <- pnorm(x)
  mod <- armacop(order = c(2, 1))
  f0 <- tsfit(u, mod, pseudo = FALSE)
  f1 <- tsfit(u, mod, init = armacop(ar = c(-0.3, 0.2), ma = -0.5), pseudo = FALSE)
  f2 <- tsfit(u, mod, init = c(0.1, 0.1, 0.1), pseudo = FALSE)
  f3 <- tsfit(u, armacop(ar = c(0.2, 0.2), ma = 0.1), pseudo = FALSE)
  expect_equal(f1@loglik, f0@loglik, tolerance = 1e-5)
  expect_equal(f2@loglik, f0@loglik, tolerance = 1e-5)
  expect_equal(f3@loglik, f0@loglik, tolerance = 1e-5)
})

test_that("white noise data fit without error and give small coefficients", {
  set.seed(10)
  u <- runif(500)
  f <- tsfit(u, armacop(order = c(1, 1)))
  expect_lt(abs(f@loglik), 5)
  f2 <- tsfit(u, armacop(order = c(1, 0)))
  expect_lt(abs(coef(f2)), 0.2)
})

test_that("NA values and ranks are handled", {
  x <- sim_arma_scores(800, ar = 0.6, seed = 11)
  x[c(10, 11, 400)] <- NA
  u <- (rank(x, na.last = "keep") - 0.5) / sum(!is.na(x))
  f <- tsfit(u, armacop(order = c(1, 0)))
  expect_equal(f@nobs, 797L)
  expect_equal(unname(coef(f)), 0.6, tolerance = 0.1)
})

test_that("input is validated", {
  mod <- armacop(order = c(1, 0))
  expect_error(tsfit("a", mod), "numeric")
  expect_error(tsfit(c(0, runif(50)), mod), "boundaryadjust")
  expect_error(tsfit(c(1, runif(50)), mod), "boundaryadjust")
  expect_error(tsfit(runif(10), mod), "at least 20")
  expect_error(tsfit(runif(100), mod, se = "hessian"), "not implemented")
  expect_error(tsfit(runif(100), mod, init = armacop(order = c(1, 1))), "unspecified|same")
  expect_error(tsfit(runif(100), mod, init = armacop(ar = 0.1, ma = 0.1)), "same orders")
  expect_error(tsfit(runif(100), mod, init = c(0.1, 0.2)), "init")
})

test_that("ts and zoo-like inputs are accepted as plain numeric", {
  set.seed(12)
  u <- ts(runif(100))
  expect_s4_class(tsfit(u, armacop(order = c(1, 0))), "armacopfit")
})

test_that("show works for fitted objects", {
  set.seed(13)
  f <- tsfit(runif(100), armacop(order = c(1, 0)))
  expect_output(show(f), "armacopfit")
  expect_output(show(f), "loglik")
})
