# ---- effective-lag rule (stubbed) -----------------------------------------

stub_rule <- function(K_seq, dll, K_start = 10L, tol = 0.01, max_refit = 3L) {
  # K_of returns the next value of K_seq at each call; ll_of returns a
  # log-likelihood that is linear in K with slope 'dll'
  calls <- 0L
  fits <- list()
  fit_fixed <- function(K, z) { fits[[length(fits) + 1L]] <<- K; list(z = z + 1) }
  K_of <- function(z) { calls <<- calls + 1L; K_seq[min(calls, length(K_seq))] }
  ll_of <- function(z, K) dll * K
  res <- dv_k_rule(fit_fixed, K_of, ll_of, K_start, z0 = 0, tol = tol, max_refit = max_refit)
  list(res = res, fits = unlist(fits))
}

test_that("K rule: unchanged K is accepted without refit", {
  r <- stub_rule(10L, 5)
  expect_equal(r$res$K, 10L)
  expect_equal(r$res$nrefit, 0L)
  expect_equal(r$fits, 10L)
})

test_that("K rule: a negligible change is accepted and the new K stored", {
  r <- stub_rule(14L, 0.002)
  expect_equal(r$res$K, 14L)
  expect_equal(r$res$nrefit, 0L)
  expect_equal(r$res$dloglik, 0.008)
})

test_that("K rule: a material increase triggers a refit with the larger K", {
  r <- stub_rule(c(20L, 20L), 0.5)
  expect_equal(r$fits, c(10L, 20L))
  expect_equal(r$res$K, 20L)
  expect_equal(r$res$nrefit, 1L)
})

test_that("K rule: K never decreases", {
  r <- stub_rule(c(20L, 12L), 0.5)
  expect_equal(r$res$K, 20L)
  expect_equal(r$fits, c(10L, 20L))
  r2 <- stub_rule(6L, 0.5) # decrease with a material change: keep old K
  expect_equal(r2$res$K, 10L)
  expect_equal(r2$res$nrefit, 0L)
  r3 <- stub_rule(6L, 0.001) # negligible: new smaller K accepted
  expect_equal(r3$res$K, 6L)
})

test_that("K rule: refits are capped with a warning", {
  expect_warning(r <- stub_rule(c(20L, 30L, 40L, 50L, 60L), 0.5, max_refit = 3L), "still changing")
  expect_equal(r$res$nrefit, 3L)
  expect_false(r$res$converged)
  expect_length(r$fits, 4L)
})

# ---- fitting ----------------------------------------------------------------

test_that("a Gaussian d-vine fit agrees with the armacop fit", {
  set.seed(1)
  u <- tssim(armacop(ar = c(0.5, -0.2)), 800)
  fg <- tsfit(u, dvinecop(armacop(order = c(2, 0))), pseudo = FALSE)
  fa <- tsfit(u, armacop(order = c(2, 0)), pseudo = FALSE)
  expect_equal(coef(fg), coef(fa), tolerance = 1e-4)
  expect_equal(fg@loglik, fa@loglik, tolerance = 1e-6)
  expect_equal(fg@model@K, 2L)
  expect_equal(fg@nrefit, 0L)
  # with an MA part: the truncation is part of the model, agreement is approximate
  u2 <- tssim(armacop(ar = 0.6, ma = 0.4), 600)
  fg2 <- tsfit(u2, dvinecop(armacop(order = c(1, 1)), tautol = 1e-6), pseudo = FALSE)
  fa2 <- tsfit(u2, armacop(order = c(1, 1)), pseudo = FALSE)
  expect_equal(unname(coef(fg2)), unname(coef(fa2)), tolerance = 5e-3)
  expect_equal(fg2@loglik, fa2@loglik, tolerance = 1e-4)
})

test_that("tsfit recovers a finite-order model with a t copula and a Gumbel copula", {
  set.seed(2)
  truth <- dvinecop(armacop(ar = c(0.5, -0.2)), pairs = list(pc("t", aux = 6), pc("gumbel")))
  u <- tssim(truth, 2000)
  f <- tsfit(u, dvinecop(armacop(order = c(2, 0)), pairs = list(pc("t", aux = 6), pc("gumbel"))),
             pseudo = FALSE)
  expect_s4_class(f, "dvinecopfit")
  expect_equal(names(coef(f)), c("ar1", "ar2", "nu1"))
  expect_equal(unname(coef(f)[1:2]), c(0.5, -0.2), tolerance = 0.1)
  expect_gt(coef(f)["nu1"], 2)
  expect_lt(coef(f)["nu1"], 20)
  expect_equal(f@npar, 3L)
  # the fit beats the starting values and a Gaussian model
  expect_gte(f@loglik, dv_model_loglik(f@start, u) - 1e-6)
  fg <- tsfit(u, dvinecop(armacop(order = c(2, 0))), pseudo = FALSE)
  expect_gt(f@loglik, fg@loglik)
  expect_lt(AIC(f), AIC(fg))
})

test_that("an ARMA model with substitutions is recovered and K is stored", {
  set.seed(3)
  specs <- list(pc("gumbel"), pc("clayton", pos = 180, neg = 270))
  truth <- dvinecop(armacop(ar = 0.6, ma = 0.4), pairs = specs)
  u <- tssim(truth, 800)
  f <- tsfit(u, dvinecop(armacop(order = c(1, 1)), pairs = specs), pseudo = FALSE)
  expect_equal(unname(coef(f)), c(0.6, 0.4), tolerance = 0.15)
  expect_gte(f@model@K, 2L)
  expect_equal(f@model@K, dv_effective_lag(
    arma_coef_expanded(f@model@arma)$ar, arma_coef_expanded(f@model@arma)$ma, 2L)$K,
    tolerance = 3)
  expect_equal(f@loglik, dv_model_loglik(f@model, u), tolerance = 1e-8)
  expect_equal(as.numeric(logLik(f)), f@loglik)
})

test_that("starting values: object, init model, numeric init", {
  set.seed(4)
  specs <- list(pc("t", aux = 6), pc("gumbel"))
  truth <- dvinecop(armacop(ar = c(0.5, -0.2)), pairs = specs)
  u <- tssim(truth, 600)
  f0 <- tsfit(u, dvinecop(armacop(order = c(2, 0)), pairs = specs), pseudo = FALSE)
  f1 <- tsfit(u, dvinecop(armacop(ar = c(0.1, 0.1)), pairs = specs), pseudo = FALSE)
  f2 <- tsfit(u, dvinecop(armacop(order = c(2, 0)), pairs = specs), pseudo = FALSE,
              init = c(0.3, 0.1, 15))
  f3 <- tsfit(u, dvinecop(armacop(order = c(2, 0)), pairs = specs), pseudo = FALSE,
              init = truth)
  for (f in list(f1, f2, f3)) expect_equal(f@loglik, f0@loglik, tolerance = 1e-3)
  # warm start from a fitted model: the model itself is a known model
  f4 <- tsfit(u, f0@model, pseudo = FALSE)
  expect_equal(f4@loglik, f0@loglik, tolerance = 1e-6)
  expect_error(tsfit(u, dvinecop(armacop(order = c(2, 0)), pairs = specs), init = c(0.3, 0.1)), "init")
  expect_error(tsfit(u, dvinecop(armacop(order = c(2, 0)), pairs = specs),
                     init = dvinecop(armacop(ar = c(0.1, 0.1)), pairs = list(pc("gumbel"), pc("gumbel")))), "init")
})

test_that("a BB1 copula starts at a feasible delta and is estimated", {
  set.seed(5)
  truth <- dvinecop(armacop(ar = 0.7), pairs = list(pc("bb1", aux = 1.4)))
  u <- tssim(truth, 800)
  f <- tsfit(u, dvinecop(armacop(order = c(1, 0)), pairs = list(pc("bb1"))), pseudo = FALSE)
  expect_equal(names(coef(f)), c("ar1", "delta1"))
  expect_gte(coef(f)["delta1"], 1)
  expect_equal(unname(coef(f)["ar1"]), 0.7, tolerance = 0.1)
  expect_error(tsfit(u, dvinecop(armacop(ar = 0.3), pairs = list(pc("bb1", aux = 3))), pseudo = FALSE),
               "infeasible")
})

test_that("tsfit validates its input", {
  set.seed(6)
  u <- tssim(armacop(ar = 0.5), 100)
  m <- dvinecop(armacop(order = c(1, 0)))
  free <- dvinecop(pairs = list(rvinecopulib::bicop_dist("gumbel", 0, 2)))
  expect_error(tsfit(u, free), "cannot be fitted")
  u2 <- u
  u2[3] <- NA
  expect_error(tsfit(u2, m), "complete")
  expect_error(tsfit(c(0, u[-1]), m), "boundaryadjust")
  expect_error(tsfit(u[1:10], m), "at least 20")
  expect_error(tsfit(u, m, se = "foo"), "se")
  expect_error(tsfit(u, m, se = "bootstrap", B = 3), "B")
})

test_that("methods on the fitted object delegate to the model", {
  set.seed(7)
  specs <- list(pc("gumbel"))
  u <- tssim(dvinecop(armacop(ar = 0.6), pairs = specs), 300)
  f <- tsfit(u, dvinecop(armacop(order = c(1, 0)), pairs = specs), pseudo = FALSE)
  expect_equal(ktau(f, 3), ktau(f@model, 3))
  expect_equal(resid(f), resid(f@model, u))
  expect_equal(residuals(f, scale = "normal"), resid(f@model, u, scale = "normal"))
  expect_equal(predict(f, type = "df", x = c(0.2, 0.6)), predict(f@model, u, type = "df", x = c(0.2, 0.6)))
  expect_length(tssim(f), 300)
  expect_length(tssim(f, 20), 20)
  expect_output(show(f), "dvinecopfit")
  ut <- ts(u, start = c(2000, 1), frequency = 12)
  ft <- tsfit(ut, dvinecop(armacop(order = c(1, 0)), pairs = specs), pseudo = FALSE)
  expect_s3_class(resid(ft), "ts")
})

test_that("pseudo-observations are accepted", {
  set.seed(8)
  specs <- list(pc("clayton", pos = 180))
  x <- qt(tssim(dvinecop(armacop(ar = 0.6), pairs = specs), 500), df = 4)
  f <- tsfit(rank(x) / 501, dvinecop(armacop(order = c(1, 0)), pairs = specs), pseudo = TRUE)
  # the pseudo-likelihood estimate is noisier than with known margins for a
  # strongly tail-dependent copula: only check that it is sensible
  expect_gt(coef(f)["ar1"], 0.3)
  expect_lt(coef(f)["ar1"], 0.9)
})

# ---- standard errors -------------------------------------------------------

test_that("Hessian standard errors are positive and carried to the coefficients", {
  set.seed(9)
  specs <- list(pc("t", aux = 6), pc("gumbel"))
  u <- tssim(dvinecop(armacop(ar = c(0.5, -0.2)), pairs = specs), 1000)
  f <- tsfit(u, dvinecop(armacop(order = c(2, 0)), pairs = specs), pseudo = FALSE, se = "hessian")
  expect_equal(f@se_method, "hessian")
  expect_true(all(f@se > 0))
  expect_equal(names(f@se), names(coef(f)))
  expect_equal(vcov(f), t(vcov(f)))
  expect_lt(f@se["ar1"], 0.1)
  expect_output(show(f), "standard errors: hessian")
})

test_that("bootstrap standard errors have the right shape and agree roughly with Hessian ones", {
  set.seed(10)
  specs <- list(pc("gumbel"))
  u <- tssim(dvinecop(armacop(ar = 0.6), pairs = specs), 600)
  h <- tsfit(u, dvinecop(armacop(order = c(1, 0)), pairs = specs), pseudo = FALSE, se = "hessian")
  b <- tsfit(u, dvinecop(armacop(order = c(1, 0)), pairs = specs), pseudo = FALSE,
             se = "bootstrap", B = 40)
  expect_equal(b@se_method, "bootstrap")
  expect_equal(dim(b@boot), c(40L, 1L))
  expect_false(anyNA(b@boot))
  expect_equal(unname(b@se), unname(h@se), tolerance = 0.5)
  expect_equal(b@estimate, h@estimate)
  r <- tsfit(u, dvinecop(armacop(order = c(1, 0)), pairs = specs), pseudo = TRUE,
             se = TRUE, B = 10)
  expect_equal(dim(r@boot), c(10L, 1L))
  r2 <- tsfit(u, dvinecop(armacop(order = c(1, 0)), pairs = specs), pseudo = FALSE,
              se = "bootstrap", B = 10, bootstart = "user")
  expect_equal(dim(r2@boot), c(10L, 1L))
  expect_error(tsfit(u, dvinecop(armacop(order = c(1, 0)), pairs = specs), se = "bootstrap",
                     B = 10, bootstart = "foo"), "arg")
})

# ---- oracle ----------------------------------------------------------------

test_that("the fit agrees with the legacy sdvinecopula fit", {
  skip_if_no_oracle()
  skip_on_cran()
  set.seed(11)
  truth <- dvinecop(armacop(ar = 0.6, ma = 0.4), pairs = list(pc("gumbel"), pc("clayton", pos = 180, neg = 270)))
  u <- tssim(truth, 400)
  f <- tsfit(u, dvinecop(armacop(order = c(1, 1)), pairs = list(pc("gumbel"), pc("clayton", pos = 180, neg = 270))),
             pseudo = FALSE)
  old <- tscopula::sdvinecopula(pars = list(ar = 0.3, ma = 0.1), family = c("gumbel", "clayton"),
                                posrot = c(0, 180), negrot = c(90, 270), tautol = 1e-4)
  fo <- tscopula::fit(old, u)
  expect_equal(unname(coef(f)), unname(coef(fo)), tolerance = 0.02)
  expect_equal(f@loglik, as.numeric(logLik(fo)), tolerance = 0.05)
})
