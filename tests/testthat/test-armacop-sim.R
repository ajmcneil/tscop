test_that("tssim returns copula data of the right length", {
  set.seed(1)
  u <- tssim(armacop(ar = 0.5, ma = 0.4), 200)
  expect_length(u, 200)
  expect_true(all(u > 0 & u < 1))
  expect_length(tssim(armacop(ar = 0.5), 1), 1)
  expect_error(tssim(armacop(ar = 0.5), 0), "positive")
  expect_error(tssim(armacop(order = c(1, 0)), 10), "unspecified")
  set.seed(2)
  a <- tssim(armacop(ar = 0.5), 20)
  set.seed(2)
  b <- tssim(armacop(ar = 0.5), 20)
  expect_equal(a, b)
})

test_that("simulation starts in the stationary distribution (no burn-in)", {
  # short paths: the joint law of (x1, x2, x3) must be the model's Toeplitz one
  check <- function(ar, ma, nrep = 20000, seed = 3) {
    set.seed(seed)
    mod <- armacop(ar = ar, ma = ma)
    X <- qnorm(t(vapply(seq_len(nrep), function(i) tssim(mod, 3), numeric(3))))
    acf <- acf_ref(ar, ma, 2)
    expect_equal(unname(cov(X)), toeplitz(acf), tolerance = 0.05, ignore_attr = TRUE)
    expect_equal(unname(colMeans(X)), rep(0, 3), tolerance = 0.05)
  }
  check(0.5, 0.4)
  check(0.95, -0.85)
  check(numeric(0), 0.6)
})

test_that("long simulated series have the model autocorrelation", {
  set.seed(4)
  mod <- armacop(ar = c(0.5, 0.2), ma = 0.4)
  x <- qnorm(tssim(mod, 2e5))
  expect_equal(var(x), 1, tolerance = 0.02)
  expect_equal(as.numeric(acf(x, lag.max = 5, plot = FALSE)$acf[-1]),
               acf_ref(c(0.5, 0.2), 0.4, 5)[-1], tolerance = 0.02)
})

test_that("seasonal models simulate with the right autocorrelation", {
  set.seed(5)
  mod <- armacop(ar = 0.3, sma = 0.5, period = 4)
  e <- arma_coef_expanded(mod)
  x <- qnorm(tssim(mod, 2e5))
  expect_equal(as.numeric(acf(x, lag.max = 8, plot = FALSE)$acf[-1]),
               acf_ref(e$ar, e$ma, 8)[-1], tolerance = 0.02)
})

test_that("tssim works for fitted objects and defaults to the data length", {
  set.seed(6)
  u <- tssim(armacop(ar = 0.6), 300)
  f <- tsfit(u, armacop(order = c(1, 0)))
  expect_length(tssim(f), 300)
  expect_length(tssim(f, 50), 50)
})

test_that("simulate-then-fit round trip", {
  set.seed(7)
  u <- tssim(armacop(ar = 0.5, ma = 0.4), 3000)
  f <- tsfit(u, armacop(order = c(1, 1)), pseudo = FALSE)
  expect_equal(unname(coef(f)), c(0.5, 0.4), tolerance = 0.1)
})

test_that("n has defaults", {
  expect_length(tssim(armacop(ar = 0.5)), 1000)
})
