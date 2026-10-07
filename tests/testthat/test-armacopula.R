test_that("known models are constructed with named coefficients", {
  m <- armacopula(ar = c(0.5, 0.2), ma = 0.4)
  expect_s4_class(m, "armacopula")
  expect_equal(m@name, "ARMA(2,1)")
  expect_equal(m@order, c(2L, 1L, 0L, 0L))
  expect_true(is.na(m@period))
  expect_equal(coef(m), c(ar1 = 0.5, ar2 = 0.2, ma1 = 0.4))
  expect_true(arma_is_specified(m))
  expect_equal(arma_coef_expanded(m), list(ar = c(0.5, 0.2), ma = 0.4))
})

test_that("seasonal models record the period and expand correctly", {
  m <- armacopula(ar = 0.5, ma = 0.4, sar = 0.3, sma = 0.2, period = 4)
  expect_equal(m@name, "SARMA(1,1)(1,1)[4]")
  expect_equal(m@period, 4)
  expect_equal(names(coef(m)), c("ar1", "ma1", "sar1", "sma1"))
  e <- arma_coef_expanded(m)
  expect_equal(e$ar, c(0.5, 0, 0, 0.3, -0.15))
  expect_equal(e$ma, c(0.4, 0, 0, 0.2, 0.08))
})

test_that("structures to be fitted hold NA coefficients", {
  s <- armacopula(order = c(1, 1))
  expect_equal(s@order, c(1L, 1L, 0L, 0L))
  expect_equal(coef(s), c(ar1 = NA_real_, ma1 = NA_real_))
  expect_false(arma_is_specified(s))
  expect_error(arma_coef_expanded(s), "unspecified")
  ss <- armacopula(order = c(1, 0), seasonal = c(0, 2), period = 12)
  expect_equal(ss@name, "SARMA(1,0)(0,2)[12]")
  expect_equal(names(coef(ss)), c("ar1", "sma1", "sma2"))
})

test_that("period is stored as NA without seasonal terms", {
  expect_true(is.na(armacopula(ar = 0.5, period = 12)@period))
  expect_true(is.na(armacopula(order = c(1, 0), seasonal = c(0, 0), period = 12)@period))
})

test_that("inadmissible and malformed input is rejected", {
  expect_error(armacopula(), "at least one")
  expect_error(armacopula(order = c(0, 0)), "at least one")
  expect_error(armacopula(ar = 1.1), "non-stationary AR")
  expect_error(armacopula(ar = 1), "non-stationary AR")
  expect_error(armacopula(ma = -1.2), "non-invertible MA")
  expect_error(armacopula(sar = 1.5, period = 4), "non-stationary seasonal AR")
  expect_error(armacopula(sma = 1.5, period = 4), "non-invertible seasonal MA")
  expect_error(armacopula(sar = 0.3), "period")
  expect_error(armacopula(sar = 0.3, period = 1), "period")
  expect_error(armacopula(sar = 0.3, period = 2.5), "period")
  expect_error(armacopula(order = c(1, 1), seasonal = c(1, 0)), "period")
  expect_error(armacopula(ar = NULL, ma = 0.5), "NULL")
  expect_error(armacopula(ar = NA_real_), "missing")
  expect_error(armacopula(ar = "a"), "numeric")
  expect_error(armacopula(ar = 0.5, order = c(1, 1)), "not both")
  expect_error(armacopula(seasonal = c(1, 0), period = 4), "order")
  expect_error(armacopula(order = c(1.5, 1)), "whole")
  expect_error(armacopula(order = c(1, -1)), "non-negative")
  expect_error(armacopula(order = 1), "two")
})

test_that("admissibility is also enforced on modified objects", {
  m <- armacopula(ar = 0.5)
  m@pars$ar <- 1.5
  expect_error(validObject(m), "non-stationary")
})

test_that("show prints the name, period and parameters", {
  expect_output(show(armacopula(ar = 0.5, ma = 0.4)), "ARMA\\(1,1\\)")
  expect_output(show(armacopula(ar = 0.5, ma = 0.4)), "ar1")
  expect_output(show(armacopula(sar = 0.3, period = 12)), "period: 12")
  expect_output(show(armacopula(order = c(1, 1))), "not yet estimated")
})
