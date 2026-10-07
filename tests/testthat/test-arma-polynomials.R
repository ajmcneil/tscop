test_that("seasonal expansion matches the hand-computed product", {
  e <- arma_expand(ar = 0.5, ma = 0.4, sar = 0.3, sma = 0.2, period = 4)
  expect_equal(e$ar, c(0.5, 0, 0, 0.3, -0.15))
  expect_equal(e$ma, c(0.4, 0, 0, 0.2, 0.08))
  e2 <- arma_expand(sar = 0.3, period = 12)
  expect_equal(e2$ar, c(rep(0, 11), 0.3))
  expect_length(e2$ma, 0)
  expect_error(arma_expand(sar = 0.3), "period")
  expect_error(arma_expand(sma = 0.3, period = 1), "period")
})

test_that("pacf maps round-trip and agree with ARMAacf", {
  ar <- c(0.5, -0.2, 0.1)
  expect_equal(ar_to_pacf(ar), as.numeric(ARMAacf(ar = ar, lag.max = 3, pacf = TRUE)))
  expect_equal(pacf_to_ar(ar_to_pacf(ar)), ar)
  set.seed(1)
  z <- rnorm(6)
  expect_equal(ar_to_free(free_to_ar(z)), z)
  ma <- c(0.3, 0.2)
  expect_equal(free_to_ma(ma_to_free(ma)), ma)
})

test_that("stationarity check agrees with polyroot", {
  set.seed(2)
  for (i in 1:200) {
    p <- sample(1:5, 1)
    ar <- runif(p, -1, 1) * runif(1, 0.3, 1.2)
    ref <- all(Mod(polyroot(c(1, -ar))) > 1)
    expect_equal(arma_is_stationary(ar), ref)
  }
  expect_true(arma_is_stationary(numeric(0)))
  expect_false(arma_is_stationary(1))
  expect_false(arma_is_invertible(-1.2))
  expect_true(arma_is_invertible(c(0.4, 0.1)))
})

test_that("any unconstrained vector gives a stationary polynomial", {
  set.seed(3)
  for (i in 1:100) {
    z <- rnorm(sample(1:8, 1))
    expect_true(arma_is_stationary(free_to_ar(z)))
    expect_true(arma_is_invertible(free_to_ma(z)))
  }
})
