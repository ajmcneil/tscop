test_that("Hessian standard errors match the AR(1) asymptotic formula", {
  set.seed(1)
  u <- tssim(armacop(ar = 0.6), 3000)
  f <- tsfit(u, armacop(order = c(1, 0)), se = "hessian", pseudo = FALSE)
  expect_equal(f@se_method, "hessian")
  phi <- unname(coef(f))
  expect_equal(unname(f@se), sqrt((1 - phi^2) / 3000), tolerance = 0.05)
  expect_equal(names(f@se), "ar1")
  expect_equal(unname(vcov(f)), matrix(f@se^2, 1, 1), ignore_attr = TRUE)
})

test_that("Hessian vcov is symmetric positive definite on a seasonal model", {
  set.seed(2)
  u <- tssim(armacop(ar = 0.4, ma = 0.3, sar = 0.5, sma = -0.3, period = 4), 1500)
  f <- tsfit(u, armacop(order = c(1, 1), seasonal = c(1, 1), period = 4),
             se = "hessian", pseudo = FALSE)
  V <- vcov(f)
  expect_equal(V, t(V))
  expect_true(all(eigen(V, symmetric = TRUE)$values > 0))
  expect_equal(dimnames(V), list(names(coef(f)), names(coef(f))))
  expect_true(all(f@se > 0 & f@se < 0.2))
})

test_that("Hessian standard errors match the empirical spread of the estimates", {
  skip_on_cran()
  set.seed(3)
  mod <- armacop(ar = 0.5, ma = 0.3)
  est <- t(vapply(1:60, function(i) {
    f <- tsfit(tssim(mod, 600), armacop(order = c(1, 1)), se = "hessian", pseudo = FALSE)
    c(coef(f), f@se)
  }, numeric(4)))
  expect_equal(unname(colMeans(est[, 3:4])), unname(apply(est[, 1:2], 2, sd)), tolerance = 0.25)
})

test_that("bootstrap standard errors agree with Hessian standard errors for PIT data", {
  set.seed(4)
  u <- tssim(armacop(ar = 0.6), 800)
  h <- tsfit(u, armacop(order = c(1, 0)), se = "hessian", pseudo = FALSE)
  b <- tsfit(u, armacop(order = c(1, 0)), se = "bootstrap", B = 80, pseudo = FALSE)
  expect_equal(b@se_method, "bootstrap")
  expect_equal(dim(b@boot), c(80L, 1L))
  expect_equal(unname(b@se), unname(h@se), tolerance = 0.3)
  expect_equal(unname(b@vcov), unname(cov(b@boot)), tolerance = 1e-8)
  expect_equal(b@estimate, h@estimate)
})

test_that("bootstrap handles ranks, missing values and TRUE", {
  set.seed(5)
  x <- qnorm(tssim(armacop(ar = 0.5), 400))
  x[c(5, 100)] <- NA
  u <- (rank(x, na.last = "keep") - 0.5) / sum(!is.na(x))
  f <- tsfit(u, armacop(order = c(1, 0)), se = TRUE, B = 20)
  expect_equal(f@se_method, "bootstrap")
  expect_equal(f@nobs, 398L)
  expect_equal(nrow(f@boot), 20L)
  expect_true(f@se > 0)
})

test_that("bootstrap reruns from the user's starting values", {
  set.seed(6)
  u <- tssim(armacop(ar = 0.5, ma = 0.3), 500)
  f <- tsfit(u, armacop(ar = 0.2, ma = 0.1), se = "bootstrap", B = 15, pseudo = FALSE)
  expect_equal(dim(f@boot), c(15L, 2L))
  expect_false(anyNA(f@boot))
})

test_that("show prints standard errors when available", {
  set.seed(7)
  u <- tssim(armacop(ar = 0.5), 300)
  f <- tsfit(u, armacop(order = c(1, 0)), se = "hessian")
  expect_output(show(f), "standard errors: hessian")
  expect_output(show(f), "se")
})
