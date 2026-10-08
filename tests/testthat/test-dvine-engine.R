test_that("engine matches the plain rvinecopulib loop for mixed families and rotations", {
  set.seed(1)
  u <- runif(60)
  cases <- list(
    list(pair("gumbel", 90, 1.8), pair("clayton", 270, 1.2), pair("joe", 180, 1.5)),
    list(pair("frank", 0, -3), pair("t", 0, c(0.4, 5)), pair("bb1", 90, c(0.4, 1.3)),
         pair("gumbel", 0, 1.2)),
    list(pair("gauss", 0, 0.5), pair("gauss", 0, -0.3), pair("gauss", 0, 0.2)),
    list(pair("gauss", 0, 0.5), pair("gumbel", 0, 1.4), pair("gauss", 0, 0.2), pair("gauss", 0, -0.1)),
    list(pair("clayton", 0, 1), pair("gauss", 0, 0.3), pair("gauss", 0, 0.2), pair("frank", 0, 1.5))
  )
  for (cs in cases) {
    expect_equal(dv_loglik(u, cs), dv_loglik_ref(u, cs), tolerance = 1e-9)
  }
})

test_that("engine matches rvinecopulib::vinecop_dist for rotated and mixed copulas", {
  set.seed(2)
  n <- 7
  cases <- list(
    list(pair("gumbel", 90, 1.8), pair("clayton", 270, 1.2)),
    list(pair("joe", 0, 2), pair("gumbel", 180, 1.5), pair("frank", 0, -3)),
    list(pair("t", 0, c(0.5, 5)), pair("gauss", 0, -0.3), pair("bb1", 90, c(0.4, 1.3))),
    list(pair("gauss", 0, 0.5), pair("frank", 0, 2), pair("gauss", 0, -0.3), pair("gumbel", 270, 1.4))
  )
  for (cs in cases) {
    for (rep in 1:3) {
      u <- runif(n)
      expect_equal(dv_loglik(u, cs), dv_loglik_vinecop(u, cs), tolerance = 1e-8)
    }
  }
})

test_that("rotations 90 and 270 are distinguished", {
  u <- pnorm(as.numeric(arima.sim(list(ar = -0.6), 100)))
  a <- dv_loglik(u, list(pair("gumbel", 90, 2)))
  b <- dv_loglik(u, list(pair("gumbel", 270, 2)))
  expect_false(isTRUE(all.equal(a, b)))
})

test_that("Gaussian pairs reproduce the exact Gaussian copula likelihood of an ARMA process", {
  set.seed(3)
  n <- 40
  for (cs in list(list(ar = 0.5, ma = 0.4), list(ar = c(0.5, 0.3), ma = numeric(0)),
                  list(ar = numeric(0), ma = 0.6), list(ar = 0.95, ma = -0.85))) {
    u <- pnorm(sim_arma_scores(n, cs$ar, cs$ma))
    alpha <- arma_pacf(cs$ar, cs$ma, n - 1)
    pairs <- lapply(alpha, function(a) pair("gauss", 0, a))
    expect_equal(dv_loglik(u, pairs), arma_copula_loglik(u, cs$ar, cs$ma), tolerance = 1e-8)
  }
})

test_that("a pure AR(p) is exactly a d-vine of order p", {
  set.seed(4)
  u <- pnorm(sim_arma_scores(100, ar = c(0.5, 0.3)))
  alpha <- arma_pacf(c(0.5, 0.3), numeric(0), 2)
  pairs <- lapply(alpha, function(a) pair("gauss", 0, a))
  expect_equal(dv_loglik(u, pairs), arma_copula_loglik(u, c(0.5, 0.3), numeric(0)), tolerance = 1e-8)
})

test_that("the engine handles short series and zero trees", {
  expect_equal(dv_loglik(c(0.3, 0.6), list()), 0)
  u <- runif(5)
  pr <- list(pair("gauss", 0, 0.4), pair("gauss", 0, 0.2), pair("gauss", 0, 0.1),
             pair("gauss", 0, 0.1), pair("gauss", 0, 0.1), pair("gauss", 0, 0.1))
  expect_equal(dv_loglik(u, pr), dv_loglik(u, pr[1:4]))
  expect_equal(dv_loglik(u, pr), dv_loglik_ref(u, pr))
})

test_that("resolved pair copulas have the requested Kendall's tau", {
  taus <- c(-0.7, -0.3, -0.05, 0, 0.05, 0.3, 0.7)
  for (tau in taus) {
    for (sp in list(pc("gauss"), pc("t", aux = 5), pc("frank"), pc("gumbel"),
                    pc("clayton", pos = 180, neg = 270), pc("joe", pos = 180))) {
      pr <- dv_resolve_pair(sp, tau)
      expect_equal(dv_pair_ktau(pr), tau, tolerance = 1e-6)
    }
  }
  # BB1 with delta = 1.5 needs |tau| >= 1 - 1/delta = 1/3
  for (tau in c(-0.7, -0.4, 0.4, 0.7)) {
    pr <- dv_resolve_pair(pc("bb1", aux = 1.5), tau)
    expect_equal(dv_pair_ktau(pr), tau, tolerance = 1e-6)
  }
  expect_null(dv_resolve_pair(pc("bb1", aux = 1.5), 0.2))
  expect_null(dv_resolve_pair(pc("bb1", aux = 1.5), -0.2))
  expect_error(dv_resolve_pair(pc("t"), 0.3), "auxiliary")
})

test_that("the rotation used follows the sign of tau", {
  sp <- pc("gumbel", pos = 180, neg = 270)
  expect_equal(dv_resolve_pair(sp, 0.3)$rotation, 180)
  expect_equal(dv_resolve_pair(sp, -0.3)$rotation, 270)
  expect_equal(dv_resolve_pair(sp, 0)$rotation, 180)
  expect_equal(dv_resolve_pair(pc("frank"), -0.3)$rotation, 0)
})

test_that("closed-form ktau_to_par agrees with rvinecopulib", {
  for (fam in c("clayton", "gumbel", "joe", "frank")) {
    for (tau in c(0.01, 0.2, 0.6)) {
      expect_equal(dv_ktau_to_par(fam, tau),
                   as.numeric(rvinecopulib::ktau_to_par(fam, tau)), tolerance = 1e-8)
    }
  }
})

# ---- effective lag -----------------------------------------------------

brute_K <- function(ar, ma, tautol, L = 3000) {
  tau <- (2 / pi) * asin(as.numeric(ARMAacf(ar, ma, lag.max = L, pacf = TRUE)))
  max(which(abs(tau) > tautol))
}

test_that("effective lag equals the brute-force last exceedance", {
  for (cs in list(list(ar = 0.5, ma = 0.4), list(ar = 0.95, ma = -0.85),
                  list(ar = numeric(0), ma = 0.8), list(ar = c(0.4, 0.2), ma = -0.5))) {
    for (tol in c(1e-3, 1e-4, 1e-6)) {
      k <- dv_effective_lag(cs$ar, cs$ma, tautol = tol)
      expect_equal(k$K, brute_K(cs$ar, cs$ma, tol))
      expect_equal(k$limit, "tol")
    }
  }
})

test_that("seasonal bumps are not cut off", {
  e <- arma_expand(ma = 0.3, sma = 0.8, period = 12)
  k <- dv_effective_lag(numeric(0), e$ma, tautol = 1e-4, period = 12)
  expect_equal(k$K, brute_K(numeric(0), e$ma, 1e-4))
  expect_gt(k$K, 24)
})

test_that("pure AR order is exact, with floor and maxlag", {
  expect_equal(dv_effective_lag(c(0.5, 0.2, 0.1), numeric(0))$K, 3L)
  expect_equal(dv_effective_lag(c(0.5, 0.2, 0.1), numeric(0))$limit, "ar")
  e <- arma_expand(ar = 0.3, sar = 0.5, period = 4)
  expect_equal(dv_effective_lag(e$ar, e$ma, period = 4)$K, 5L)
  k <- dv_effective_lag(0.5, numeric(0), K0 = 3)
  expect_equal(c(k$K, k$limit), c(3L, "floor"), ignore_attr = TRUE)
  k <- dv_effective_lag(0.5, 0.4, maxlag = 4)
  expect_equal(k$K, 4L)
  expect_equal(k$limit, "maxlag")
  k <- dv_effective_lag(0.5, 0.4, K0 = 6, maxlag = 4)
  expect_equal(k$K, 6L)
  expect_equal(k$limit, "floor")
})

test_that("a kpacf that does not decay triggers the cap warning", {
  expect_warning(k <- dv_effective_lag(0.5, 0.999, cap = 100), "still above")
  expect_equal(k$limit, "cap")
})

test_that("arma_pacf is accurate far out for a slowly decaying model", {
  a <- arma_pacf(0.98, -0.5, 400)
  ref <- as.numeric(ARMAacf(0.98, -0.5, lag.max = 400, pacf = TRUE))
  expect_equal(a, ref, tolerance = 1e-8)
})
