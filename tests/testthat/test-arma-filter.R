# stats::ARMAacf rejects the empty model; white noise has acf (1, 0, 0, ...).
acf_ref <- function(ar, ma, lag.max) {
  if (length(ar) == 0 && length(ma) == 0) return(c(1, numeric(lag.max)))
  as.numeric(ARMAacf(ar = ar, ma = ma, lag.max = lag.max))
}

# Exact Gaussian log-likelihood from the Toeplitz correlation matrix.
toeplitz_loglik <- function(x, ar, ma) {
  ok <- !is.na(x)
  n <- length(x)
  acf <- acf_ref(ar, ma, n - 1)
  R <- toeplitz(acf)[ok, ok, drop = FALSE]
  L <- chol(R)
  z <- backsolve(L, x[ok], transpose = TRUE)
  -0.5 * (sum(ok) * log(2 * pi) + 2 * sum(log(diag(L))) + sum(z^2))
}

models <- list(
  ar1 = list(ar = 0.6, ma = numeric(0)),
  ar2 = list(ar = c(0.5, 0.3), ma = numeric(0)),
  ma1 = list(ar = numeric(0), ma = 0.5),
  ma2 = list(ar = numeric(0), ma = c(0.4, -0.3)),
  arma11 = list(ar = 0.5, ma = 0.4),
  arma21 = list(ar = c(0.5, 0.4), ma = -0.8),
  nearunit = list(ar = 0.95, ma = -0.85),
  wn = list(ar = numeric(0), ma = numeric(0))
)

test_that("autocovariances match ARMAacf and closed forms", {
  for (mod in models) {
    g <- arma_acvf(mod$ar, mod$ma, lag.max = 15)
    expect_equal(g / g[1], acf_ref(mod$ar, mod$ma, 15))
  }
  phi <- 0.5
  th <- 0.4
  expect_equal(arma_acvf(phi, th, 0), (1 + 2 * phi * th + th^2) / (1 - phi^2))
  expect_equal(arma_acvf(numeric(0), th, 1), c(1 + th^2, th))
  expect_equal(arma_acvf(phi, numeric(0), 3), phi^(0:3) / (1 - phi^2))
})

test_that("the stationary covariance agrees with the Kronecker solution", {
  mod <- models$arma21
  m <- 2
  Tm <- cbind(c(mod$ar), rbind(1, 0))
  R <- c(1, mod$ma)
  RR <- tcrossprod(R)
  kron <- matrix(solve(diag(m^2) - kronecker(Tm, Tm), as.vector(RR)), m, m)
  expect_equal(arma_stationary_cov(Tm, RR), kron)
})

test_that("filter log-likelihood matches the exact Toeplitz likelihood", {
  set.seed(10)
  x <- rnorm(80)
  for (mod in models) {
    expect_equal(arma_filter(x, mod$ar, mod$ma)$loglik,
                 toeplitz_loglik(x, mod$ar, mod$ma), tolerance = 1e-9)
  }
})

test_that("the first predictive variance is one (unit-variance constraint)", {
  for (mod in models) {
    f <- arma_filter(rnorm(20), mod$ar, mod$ma)
    expect_equal(f$F[1], 1)
    expect_equal(f$mu[1], 0)
  }
})

test_that("filter log-likelihood matches stats::arima at the constrained variance", {
  set.seed(11)
  x <- rnorm(120)
  n <- length(x)
  for (nm in c("ar2", "ma2", "arma11", "arma21")) {
    mod <- models[[nm]]
    p <- length(mod$ar)
    q <- length(mod$ma)
    fit <- arima(x, order = c(p, 0, q), include.mean = FALSE,
                 fixed = c(mod$ar, mod$ma), transform.pars = FALSE, method = "ML")
    # recover sum(log f_t) and S from arima's concentrated likelihood, then
    # evaluate at sigma^2 = 1 / gamma_0
    s2 <- fit$sigma2
    sumlogf <- -2 * fit$loglik - n * log(2 * pi * s2) - n
    g0 <- arma_acvf(mod$ar, mod$ma, 0)
    ll <- -0.5 * n * log(2 * pi / g0) - 0.5 * sumlogf - 0.5 * n * s2 * g0
    expect_equal(arma_filter(x, mod$ar, mod$ma)$loglik, ll, tolerance = 1e-6)
  }
})

test_that("missing values give the marginal Gaussian likelihood", {
  set.seed(12)
  x <- rnorm(60)
  x[c(5, 6, 30, 31, 32, 60)] <- NA
  for (mod in models) {
    expect_equal(arma_filter(x, mod$ar, mod$ma)$loglik,
                 toeplitz_loglik(x, mod$ar, mod$ma), tolerance = 1e-9)
  }
})

test_that("steady-state switching does not change the likelihood", {
  set.seed(13)
  x <- rnorm(300)
  for (mod in models) {
    a <- arma_filter(x, mod$ar, mod$ma, tol = 1e-12)$loglik
    b <- arma_filter(x, mod$ar, mod$ma, tol = -1)$loglik  # never switches
    expect_equal(a, b, tolerance = 1e-9)
  }
})

test_that("seasonal models agree with the expanded ARMA and Toeplitz likelihood", {
  set.seed(14)
  x <- rnorm(100)
  e <- arma_expand(ar = 0.4, ma = 0.3, sar = 0.5, sma = -0.4, period = 12)
  expect_equal(arma_filter(x, e$ar, e$ma)$loglik,
               toeplitz_loglik(x, e$ar, e$ma), tolerance = 1e-8)
})

test_that("copula log-likelihood subtracts the normal score density", {
  set.seed(15)
  u <- runif(50)
  x <- qnorm(u)
  expect_equal(arma_copula_loglik(u, 0.5, numeric(0)),
               toeplitz_loglik(x, 0.5, numeric(0)) - sum(dnorm(x, log = TRUE)),
               tolerance = 1e-9)
  # independence: copula density is one
  expect_equal(arma_copula_loglik(u), 0, tolerance = 1e-9)
})
