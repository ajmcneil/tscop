# exact PACF by solving the Toeplitz system for each order (independent of the
# Durbin-Levinson recursion in the package)
pacf_by_solve <- function(rho, K) {
  vapply(seq_len(K), function(k) {
    R <- toeplitz(c(1, rho[seq_len(k - 1L)]))
    phi <- solve(R, rho[seq_len(k)])
    phi[k]
  }, numeric(1))
}

test_that("constructors build known models and structures", {
  g <- arfactcop(lambda = c(0.7, 0.98), w = c(0.4, 0.3))   # unsorted input
  expect_s4_class(g, "arfactcop")
  expect_s4_class(g, "gausscop")
  expect_equal(g@lambda, c(0.98, 0.7))
  expect_equal(g@w, c(0.3, 0.4))
  expect_equal(g@name, "ARFACT(2)")
  expect_equal(names(coef(g)), c("lambda1", "lambda2", "w1", "w2"))
  expect_equal(unname(coef(g)), c(0.98, 0.7, 0.3, 0.4))
  s <- arfactcop(order = 2)
  expect_false(gc_is_specified(s))
  expect_true(gc_is_specified(g))
  nn <- arfactcop(lambda = c(0.9, 0.5), w = c(0.6, 0.4), noise = FALSE)
  expect_equal(names(coef(nn)), c("lambda1", "lambda2", "w1"))
  expect_equal(nn@name, "ARFACT(2) without noise")
  expect_output(show(g), "noise share w0: 0.3")
  expect_output(show(s), "not yet estimated")
})

test_that("invalid models are rejected", {
  expect_error(arfactcop(), "specify")
  expect_error(arfactcop(lambda = 1.2, w = 0.5), "lie in")
  expect_error(arfactcop(lambda = 0, w = 0.5), "lie in")
  expect_error(arfactcop(lambda = 0.5, w = -0.1), "non-negative")
  expect_error(arfactcop(lambda = c(0.5, 0.4), w = c(0.7, 0.6)), "at most one")
  expect_error(arfactcop(lambda = c(0.5, 0.4), w = 0.5), "same length")
  expect_error(arfactcop(lambda = c(0.5, 0.5), w = c(0.2, 0.2)), "strictly decreasing")
  expect_error(arfactcop(lambda = 0.5, w = 0.5, noise = FALSE), "sum to one")
  expect_error(arfactcop(lambda = 0.5, w = 0.5, order = 1), "not both")
  expect_error(arfactcop(order = 0), "positive")
  expect_error(arfactcop(lambda = NA_real_, w = 0.5), "missing")
})

test_that("autocorrelations and the PACF agree with an independent computation", {
  g <- arfactcop(lambda = c(0.98, 0.7), w = c(0.3, 0.4))
  rho <- 0.3 * 0.98^(1:12) + 0.4 * 0.7^(1:12)
  expect_equal(afc_acf(g@lambda, g@w, 12), rho)
  expect_equal(gc_pacf(g, 12), pacf_by_solve(rho, 12), tolerance = 1e-8)
  expect_equal(gc_tau(g, 5), (2 / pi) * asin(pacf_by_solve(rho, 5)), tolerance = 1e-8)
  # PACF values of the worked example in the theory note
  expect_equal(round(gc_pacf(g, 8), 3), c(0.574, 0.231, 0.116, 0.072, 0.055, 0.046, 0.041, 0.037))
})

test_that("an AR(1) without noise has PACF (lambda, 0, 0, ...) and one factor plus noise is an ARMA(1,1)", {
  ar1 <- arfactcop(lambda = 0.8, w = 1, noise = FALSE)
  expect_equal(gc_pacf(ar1, 5), c(0.8, 0, 0, 0, 0), tolerance = 1e-12)
  expect_equal(gc_effective_lag(ar1, 0L, 1e-4, Inf), list(K = 1L, limit = "ar"))
  lam <- 0.9
  w <- 0.6
  b <- (2 * w * lam^2 - 1 - lam^2) / (lam * (w - 1))
  r <- Re(polyroot(c(1, b, 1)))
  theta <- r[abs(r) < 1]
  expect_equal(gc_pacf(arfactcop(lambda = lam, w = w), 20), arma_pacf(lam, theta, 20), tolerance = 1e-8)
  # the closed-form inverse map
  w_back <- (1 + lam * theta) * (lam + theta) / (lam * (1 + 2 * lam * theta + theta^2))
  expect_equal(w_back, w, tolerance = 1e-10)
})

test_that("the partial autocorrelations are non-negative in random draws", {
  set.seed(1)
  worst <- Inf
  for (i in 1:1500) {
    m <- sample(1:4, 1)
    noise <- runif(1) < 0.5
    lam <- sort(runif(m, 0, 0.995), decreasing = TRUE)
    w <- rexp(m)
    w <- w / sum(w) * (if (noise) runif(1, 0.2, 1) else 1)
    g <- arfactcop(lambda = lam, w = w, noise = noise)
    worst <- min(worst, min(gc_pacf(g, 80)))
  }
  expect_gte(worst, -1e-9)
})

test_that("free-scale maps are consistent and invertible", {
  for (cs in list(list(m = 1L, noise = TRUE), list(m = 1L, noise = FALSE),
                  list(m = 2L, noise = TRUE), list(m = 3L, noise = FALSE))) {
    G <- arfactcop(order = cs$m, noise = cs$noise)
    ops <- gc_free_ops(G)
    expect_equal(ops$nfree, if (cs$noise) 2 * cs$m else 2 * cs$m - 1)
    set.seed(2)
    z <- rnorm(ops$nfree, sd = 0.8)
    model <- ops$build(z)
    expect_true(gc_is_specified(model))
    expect_equal(ops$to_free(model), z, tolerance = 1e-8)
    expect_equal(unname(ops$coef(z)), unname(coef(model)), tolerance = 1e-10)
    expect_equal(ops$pacf(z, 10), gc_pacf(model, 10), tolerance = 1e-10)
    expect_equal(ops$effective_lag(z, 1L, 1e-4, Inf), gc_effective_lag(model, 1L, 1e-4, Inf))
    # rates decreasing and weights on the simplex
    p <- afc_free_to_par(z, cs$m, cs$noise)
    expect_true(all(diff(p$lambda) < 0))
    expect_equal(sum(p$w) + p$w0, 1, tolerance = 1e-12)
    # building from natural coefficients reproduces the model
    back <- gc_from_coef(G, coef(model))
    expect_equal(unname(coef(back)), unname(coef(model)), tolerance = 1e-10)
  }
  expect_equal(gc_structure(arfactcop(order = 2)), gc_structure(arfactcop(lambda = c(.9, .5), w = c(.2, .2))))
  expect_false(identical(gc_structure(arfactcop(order = 2)), gc_structure(arfactcop(order = 2, noise = FALSE))))
})

test_that("the effective lag equals the brute-force last exceedance", {
  for (g in list(arfactcop(lambda = c(0.98, 0.7), w = c(0.3, 0.4)),
                 arfactcop(lambda = 0.9, w = 0.6),
                 arfactcop(lambda = c(0.95, 0.5, 0.2), w = c(0.2, 0.3, 0.2), noise = TRUE))) {
    for (tol in c(1e-3, 1e-4)) {
      k <- gc_effective_lag(g, 0L, tol, Inf)
      tau <- (2 / pi) * asin(gc_pacf(g, 1500))
      expect_equal(k$K, max(which(abs(tau) > tol)))
      expect_equal(k$limit, "tol")
    }
  }
  g <- arfactcop(lambda = 0.9, w = 0.6)
  expect_equal(gc_effective_lag(g, 0L, 1e-4, 5), list(K = 5L, limit = "maxlag"))
  expect_equal(gc_effective_lag(g, 12L, 1e-4, 5), list(K = 12L, limit = "floor"))
})

# ---- Kalman filter ---------------------------------------------------------

toeplitz_ll <- function(x, rho) {
  ok <- !is.na(x)
  R <- toeplitz(c(1, rho[seq_len(length(x) - 1L)]))[ok, ok, drop = FALSE]
  L <- chol(R)
  z <- backsolve(L, x[ok], transpose = TRUE)
  -0.5 * (sum(ok) * log(2 * pi) + 2 * sum(log(diag(L))) + sum(z^2))
}

test_that("the Kalman filter gives the exact Gaussian likelihood", {
  set.seed(3)
  x <- rnorm(60)
  for (cs in list(list(lambda = 0.9, w = 0.6), list(lambda = c(0.97, 0.6), w = c(0.3, 0.3)),
                  list(lambda = c(0.9, 0.5, 0.2), w = c(0.2, 0.3, 0.1)),
                  list(lambda = c(0.9, 0.4), w = c(0.5, 0.5)))) {
    w0 <- 1 - sum(cs$w)
    rho <- afc_acf(cs$lambda, cs$w, 59)
    f <- afc_filter(x, cs$lambda, cs$w, w0)
    expect_equal(f$loglik, toeplitz_ll(x, rho), tolerance = 1e-9)
    expect_equal(f$F[1], 1, tolerance = 1e-12)
    xm <- x
    xm[c(4, 5, 30, 60)] <- NA
    expect_equal(afc_filter(xm, cs$lambda, cs$w, w0)$loglik, toeplitz_ll(xm, rho), tolerance = 1e-9)
    expect_equal(afc_filter(x, cs$lambda, cs$w, w0, tol = -1)$loglik, f$loglik, tolerance = 1e-9)
  }
})

test_that("one factor plus noise has the same likelihood as the equivalent ARMA(1,1)", {
  set.seed(4)
  lam <- 0.9
  w <- 0.6
  b <- (2 * w * lam^2 - 1 - lam^2) / (lam * (w - 1))
  r <- Re(polyroot(c(1, b, 1)))
  theta <- r[abs(r) < 1]
  u <- pnorm(sim_arma_scores(200, lam, theta))
  expect_equal(afc_copula_loglik(u, lam, w, 1 - w), arma_copula_loglik(u, lam, theta), tolerance = 1e-8)
})

# ---- maximum likelihood ----------------------------------------------------

test_that("the Gaussian MLE recovers a factor model and a poor start converges to the same optimum", {
  x <- sim_arfact_scores(2500, c(0.97, 0.6), c(0.3, 0.3), seed = 5)
  u <- pnorm(x)
  fit <- arfactcop_mle(u, arfactcop(order = 2), NULL, list())
  expect_true(gc_is_specified(fit$model))
  expect_equal(unname(coef(fit$model)[1:2]), c(0.97, 0.6), tolerance = 0.15)
  fit2 <- arfactcop_mle(u, arfactcop(order = 2), arfactcop(lambda = c(0.5, 0.2), w = c(0.1, 0.1)), list())
  expect_equal(fit2$opt$value, fit$opt$value, tolerance = 1e-5)
  # one factor
  x1 <- sim_arfact_scores(2000, 0.9, 0.6, seed = 6)
  f1 <- arfactcop_mle(pnorm(x1), arfactcop(order = 1), NULL, list())
  expect_equal(unname(coef(f1$model)), c(0.9, 0.6), tolerance = 0.1)
  # no noise
  x2 <- sim_arfact_scores(1500, 0.8, 1, seed = 7)
  f2 <- arfactcop_mle(pnorm(x2), arfactcop(order = 1, noise = FALSE), NULL, list())
  expect_equal(unname(coef(f2$model)), 0.8, tolerance = 0.08)
  # prefit through the interface
  g <- gc_prefit(arfactcop(order = 1), pnorm(x1), list())
  expect_true(gc_is_specified(g))
})

# ---- use as a kpacf source in sdvinecop -----------------------------------------

test_that("sdvinecop accepts an arfactcop and resolves pair copulas from its kpacf", {
  g <- arfactcop(lambda = c(0.98, 0.7), w = c(0.3, 0.4))
  m <- sdvinecop(g, pairs = list(pc("ast"), pc("gumbel")), tail = "ast")
  expect_equal(m@mode, "tied")
  expect_gt(m@K, 5L)
  expect_output(show(m), "kpacf: 2-factor AR\\(1\\) model")
  tau <- dv_kpacf_tau(g, m@K)
  pr <- dv_resolve(m)
  for (j in seq_len(m@K)) expect_lt(abs(dv_pair_ktau(pr[[j]]) - tau[j]), 1e-6)
  expect_equal(names(coef(m)), c("lambda1", "lambda2", "w1", "w2"))
  # the structure only gives an undetermined effective lag
  expect_true(is.na(sdvinecop(arfactcop(order = 2), pairs = list(pc("ast")))@K))
})

test_that("a d-vine with an arfactcop kpacf equals the ARMA-based one for one factor plus noise", {
  set.seed(8)
  lam <- 0.9
  w <- 0.6
  b <- (2 * w * lam^2 - 1 - lam^2) / (lam * (w - 1))
  r <- Re(polyroot(c(1, b, 1)))
  theta <- r[abs(r) < 1]
  u <- tssim(armacop(ar = lam, ma = theta), 300)
  ma <- sdvinecop(armacop(ar = lam, ma = theta), pairs = list(pc("gumbel"), pc("ast")), tautol = 1e-8)
  mf <- sdvinecop(arfactcop(lambda = lam, w = w), pairs = list(pc("gumbel"), pc("ast")), tautol = 1e-8)
  expect_equal(dv_model_loglik(mf, u), dv_model_loglik(ma, u), tolerance = 1e-6)
  # Gaussian d-vine log-likelihood approaches the exact factor likelihood
  g <- sdvinecop(arfactcop(lambda = lam, w = w), tautol = 1e-12)
  expect_equal(dv_model_loglik(g, u), afc_copula_loglik(u, lam, w, 1 - w), tolerance = 1e-6)
})

test_that("free mode accepts an arfactcop tail", {
  g <- arfactcop(lambda = 0.9, w = 0.6)
  b <- list(rvinecopulib::bicop_dist("gumbel", 0, 2), rvinecopulib::bicop_dist("clayton", 0, 1))
  m <- sdvinecop(g, pairs = b)
  expect_equal(m@mode, "free")
  expect_gt(m@K, 2L)
  pr <- dv_resolve(m)
  alpha <- gc_pacf(g, m@K)
  expect_equal(pr[[4]]$parameters, alpha[4])
  u <- tssim(m, 100)
  expect_true(all(u > 0 & u < 1))
})

test_that("tsfit and autodvinecop work with an arfactcop kpacf", {
  set.seed(9)
  truth <- sdvinecop(arfactcop(lambda = 0.9, w = 0.6), pairs = list(pc("gumbel")), tautol = 1e-3)
  u <- tssim(truth, 700)
  spec <- sdvinecop(arfactcop(order = 1), pairs = list(pc("gumbel")), tautol = 1e-3)
  f <- tsfit(u, spec, pseudo = FALSE)
  expect_s4_class(f, "sdvinecopfit")
  expect_equal(names(coef(f)), c("lambda1", "w1"))
  expect_equal(unname(coef(f)), c(0.9, 0.6), tolerance = 0.2)
  expect_equal(f@npar, 2L)
  gfit <- tsfit(u, sdvinecop(arfactcop(order = 1), tautol = 1e-3), pseudo = FALSE)
  expect_gt(f@loglik, gfit@loglik)
  a <- autodvinecop(u, arfactcop(order = 1), nreplace = 2, pseudo = FALSE, verbose = FALSE, tautol = 1e-3)
  expect_s4_class(a, "sdvinecopfit")
  expect_equal(a@path$choice[1], "gumbel0")
  # standard errors
  fh <- tsfit(u, spec, pseudo = FALSE, se = "hessian")
  expect_true(all(fh@se > 0))
  # resid, predict and ktau on the fitted model
  expect_length(resid(f), 700)
  expect_length(ktau(f, 5), 5)
  expect_length(predict(f, type = "df", x = c(0.2, 0.7)), 2)
})

# ---- conversion to ARMA -----------------------------------------------------

test_that("arfact2arma reproduces the autocorrelations, PACF and likelihood", {
  set.seed(11)
  for (i in 1:25) {
    m <- sample(1:3, 1)
    noise <- runif(1) < 0.6
    lam <- sort(runif(m, 0.1, 0.99), decreasing = TRUE)
    w <- rexp(m)
    w <- w / sum(w) * (if (noise) runif(1, 0.3, 0.9) else 1)
    g <- arfactcop(lambda = lam, w = w, noise = noise)
    a <- arfact2arma(g)
    expect_s4_class(a, "armacop")
    e <- arma_coef_expanded(a)
    expect_equal(length(e$ar), m)
    expect_equal(length(e$ma), if (noise) m else m - 1L)
    expect_equal(as.numeric(ARMAacf(e$ar, e$ma, 30))[-1], afc_acf(lam, w, 30), tolerance = 1e-7)
    expect_equal(gc_pacf(a, 15), gc_pacf(g, 15), tolerance = 1e-7)
  }
  g <- arfactcop(lambda = c(0.9, 0.5), w = c(0.3, 0.4))
  u <- tssim(g, 100)
  expect_equal(arma_copula_loglik(u, arma_coef_expanded(arfact2arma(g))$ar, arma_coef_expanded(arfact2arma(g))$ma),
               afc_copula_loglik(u, g@lambda, g@w, 1 - sum(g@w)), tolerance = 1e-6)
})

test_that("arfact2arma matches the closed forms and the worked examples", {
  # one factor plus noise: ARMA(1,1) with the closed-form MA coefficient
  lam <- 0.9
  w <- 0.6
  b <- (2 * w * lam^2 - 1 - lam^2) / (lam * (w - 1))
  r <- Re(polyroot(c(1, b, 1)))
  theta <- r[abs(r) < 1]
  a <- arfact2arma(arfactcop(lambda = lam, w = w))
  expect_equal(unname(coef(a)), c(lam, theta), tolerance = 1e-8)
  # a GARCH(1,1)-type model: alpha1 = 0.1, beta1 = 0.85
  al <- 0.1
  be <- 0.85
  lam2 <- al + be
  w2 <- al * (1 - be * lam2) / (lam2 * (1 - 2 * be * lam2 + be^2))
  a2 <- arfact2arma(arfactcop(lambda = lam2, w = w2))
  expect_equal(unname(coef(a2)), c(lam2, -be), tolerance = 1e-7)
  # the two-factor example of the theory note
  a3 <- arfact2arma(arfactcop(lambda = c(0.98, 0.7), w = c(0.3, 0.4)))
  expect_equal(unname(coef(a3)), c(1.68, -0.686, -1.29224, 0.33577), tolerance = 1e-3)
  # no noise, one factor: an AR(1)
  a4 <- arfact2arma(arfactcop(lambda = 0.8, w = 1, noise = FALSE))
  expect_equal(unname(coef(a4)), 0.8)
  # coercion and fitted objects
  g <- arfactcop(lambda = c(0.98, 0.7), w = c(0.3, 0.4))
  expect_equal(coef(as(g, "armacop")), coef(arfact2arma(g)))
  f <- tsfit(tssim(g, 600), arfactcop(order = 2), pseudo = FALSE)
  expect_equal(coef(arfact2arma(f)), coef(arfact2arma(f@model)))
  expect_error(arfact2arma(arfactcop(order = 2)), "unspecified")
  expect_error(arfact2arma(1), "arfactcop")
})
