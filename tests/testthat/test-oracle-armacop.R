# Comparisons with the legacy tscopula package (numerical oracle only).

test_that("log-likelihood agrees with tscopula at fixed parameters", {
  skip_if_no_oracle()
  x <- sim_arma_scores(300, ar = c(0.5, 0.3), ma = 0.4, seed = 21)
  u <- pnorm(x)
  old <- tscopula:::armacopula_objective(c(0.5, 0.3, 0.4), c(2, 1), u)
  expect_equal(-arma_copula_loglik(u, c(0.5, 0.3), 0.4), old, tolerance = 1e-8)
  # seasonal
  e <- arma_expand(ar = 0.4, ma = 0.3, sar = 0.5, sma = -0.3, period = 4)
  olds <- tscopula:::sarmacopula_objective(c(0.4, 0.3, 0.5, -0.3), c(1, 1, 1, 1, 4), u)
  expect_equal(-arma_copula_loglik(u, e$ar, e$ma), olds, tolerance = 1e-8)
})

test_that("fitted estimates agree with tscopula", {
  skip_if_no_oracle()
  x <- sim_arma_scores(1000, ar = 0.5, ma = 0.4, seed = 22)
  u <- pnorm(x)
  new <- tsfit(u, armacop(order = c(1, 1)), pseudo = FALSE)
  old <- tscopula::fit(tscopula::armacopula(list(ar = 0.3, ma = 0.1)), u)
  expect_equal(unname(coef(new)), unname(coef(old)), tolerance = 1e-3)
  expect_equal(new@loglik, as.numeric(logLik(old)), tolerance = 1e-6)
})
