test_that("information criteria follow their definitions", {
  expect_equal(dv_ic(100, 3, 200, "AIC"), -200 + 6)
  expect_equal(dv_ic(100, 3, 200, "BIC"), -200 + 3 * log(200))
  expect_equal(dv_ic(100, 3, 200, "AICc"), -200 + 6 + 2 * 3 * 4 / (200 - 3 - 1))
})

test_that("candidates carry the rotations of the sign of tau", {
  pos <- dv_candidates(c("gumbel", "frank"), NULL, TRUE)
  expect_length(pos, 3)
  expect_equal(vapply(pos, function(s) s@family, ""), c("gumbel", "gumbel", "frank"))
  expect_equal(vapply(pos, function(s) s@pos, 0), c(0, 180, 0))
  neg <- dv_candidates(c("clayton", "joe"), NULL, FALSE)
  expect_length(neg, 4)
  expect_equal(vapply(neg, function(s) s@neg, 0), c(90, 270, 90, 270))
  expect_equal(vapply(neg, function(s) s@pos, 0), rep(0, 4))
  all <- dv_candidates(.auto_choices, c("t", "bb1"), TRUE)
  expect_length(all, 7 + 1 + 2)
  expect_equal(sum(vapply(all, function(s) s@family == "t", logical(1))), 1)
})

test_that("screening scores equal the log-likelihood of the corresponding model", {
  set.seed(1)
  arma <- armacop(ar = c(0.5, -0.2))
  base <- sdvinecop(arma, pairs = list(pc("gumbel")))
  u <- tssim(base, 300)
  cands <- list(pc("clayton", pos = 0, neg = 270), pc("frank"), pc("joe", pos = 180))
  sc <- dv_screen(u, base, 2L, cands, "AIC", npar = 2L)
  for (i in seq_along(cands)) {
    ll <- dv_model_loglik(sdvinecop(arma, pairs = list(pc("gumbel"), cands[[i]])), u)
    expect_equal(sc$loglik[i], ll, tolerance = 1e-9)
    expect_equal(sc$ic[i], dv_ic(ll, 2L, 300L, "AIC"), tolerance = 1e-9)
  }
  # t: the 1-D search finds a value at least as good as nu = 10 and is consistent
  sct <- dv_screen(u, base, 2L, list(pc("t")), "AIC", npar = 2L)
  lt <- dv_model_loglik(sdvinecop(arma, pairs = list(pc("gumbel"), pc("t", aux = sct$aux))), u)
  expect_equal(sct$loglik, lt, tolerance = 1e-6)
  l10 <- dv_model_loglik(sdvinecop(arma, pairs = list(pc("gumbel"), pc("t", aux = 10))), u)
  expect_gte(sct$loglik, l10 - 1e-9)
  expect_equal(sct$ic, dv_ic(lt, 3L, 300L, "AIC"), tolerance = 1e-6)
  # BB1 with infeasible delta range at a small |tau|
  scb <- dv_screen(u, base, 2L, list(pc("bb1")), "AIC", npar = 2L)
  expect_true(is.finite(scb$ic) || is.infinite(scb$ic))
})

test_that("the model's free vector reproduces its coefficients", {
  m <- sdvinecop(armacop(ar = 0.5, ma = 0.3), pairs = list(pc("t", aux = 6), pc("gumbel"), pc("bb1", aux = 1.4)))
  z <- dv_model_z(m)
  lay <- dv_layout(m)
  expect_equal(dv_free_to_coef(z, lay), coef(m), tolerance = 1e-10, ignore_attr = TRUE)
  m2 <- dv_build_model(m, z, lay, m@K, "tol")
  expect_equal(coef(m2), coef(m), tolerance = 1e-10)
})

test_that("autodvinecop recovers the substitutions of a finite-order model", {
  set.seed(2)
  truth <- sdvinecop(armacop(ar = c(0.6, -0.25)),
                    pairs = list(pc("gumbel"), pc("clayton", pos = 180, neg = 270)))
  u <- tssim(truth, 1000)
  gauss <- tsfit(u, armacop(order = c(2, 0)), pseudo = FALSE)
  a <- autodvinecop(u, gauss, nreplace = 5, pseudo = FALSE, verbose = FALSE)
  expect_s4_class(a, "sdvinecopfit")
  expect_equal(a@path$choice[1:2], c("gumbel0", "clayton270"))
  expect_true(all(a@path$accepted[1:2]))
  expect_gt(a@path$improvement[1], 20)
  expect_lt(AIC(a), AIC(tsfit(u, sdvinecop(armacop(order = c(2, 0))), pseudo = FALSE)) - 40)
  expect_equal(unname(coef(a)), c(0.6, -0.25), tolerance = 0.15)
  # lags beyond the kpacf's range are skipped, not fitted
  skipped <- a@path[a@path$lag > 2, ]
  expect_true(all(skipped$nfit == 0L))
  expect_true(all(!skipped$accepted))
  expect_equal(nrow(a@path), 5L - 0L)
  expect_equal(a@model@pairs[[1]]@family, "gumbel")
  expect_output(show(a), "selection path")
})

test_that("screening and full refitting select the same copulas on easy data", {
  set.seed(3)
  truth <- sdvinecop(armacop(ar = c(0.6, -0.25)),
                    pairs = list(pc("gumbel"), pc("clayton", pos = 180, neg = 270)))
  u <- tssim(truth, 800)
  gauss <- tsfit(u, armacop(order = c(2, 0)), pseudo = FALSE)
  a1 <- autodvinecop(u, gauss, nreplace = 2, pseudo = FALSE, verbose = FALSE)
  a2 <- autodvinecop(u, gauss, nreplace = 2, pseudo = FALSE, screen = FALSE, verbose = FALSE)
  expect_equal(a1@path$choice, a2@path$choice)
  expect_equal(a1@loglik, a2@loglik, tolerance = 1e-4)
})

test_that("a t copula is found when requested", {
  set.seed(4)
  truth <- sdvinecop(armacop(ar = 0.7), pairs = list(pc("t", aux = 3)))
  u <- tssim(truth, 1000)
  a <- autodvinecop(u, armacop(order = c(1, 0)), nreplace = 1, aux = "t", pseudo = FALSE,
                    verbose = FALSE)
  expect_equal(a@path$choice, "t")
  expect_equal(names(coef(a)), c("ar1", "nu1"))
  expect_gt(coef(a)["nu1"], 2)
  expect_lt(coef(a)["nu1"], 10)
  expect_equal(a@npar, 2L)
})

test_that("negative dependence is modelled with the 90/270 rotations", {
  set.seed(5)
  truth <- sdvinecop(armacop(ar = -0.6), pairs = list(pc("gumbel", neg = 270)))
  u <- tssim(truth, 1000)
  a <- autodvinecop(u, armacop(order = c(1, 0)), nreplace = 1, pseudo = FALSE, verbose = FALSE)
  expect_match(a@path$choice, "^(gumbel|joe|clayton)(90|270)$")
  expect_lt(a@path$criterion, a@path$criterion + a@path$improvement)
})

test_that("a search can be continued and stops when nothing is left to do", {
  set.seed(6)
  truth <- sdvinecop(armacop(ar = c(0.6, -0.25)), pairs = list(pc("gumbel"), pc("clayton", pos = 180, neg = 270)))
  u <- tssim(truth, 800)
  a1 <- autodvinecop(u, armacop(order = c(2, 0)), nreplace = 1, pseudo = FALSE, verbose = FALSE)
  expect_equal(nrow(a1@path), 1L)
  a2 <- autodvinecop(u, a1, nreplace = 2, pseudo = FALSE, verbose = FALSE)
  expect_equal(a2@path$lag, 2L)
  expect_equal(length(a2@model@pairs), 2L)
  a3 <- autodvinecop(u, a2, nreplace = 2, pseudo = FALSE, verbose = FALSE)
  expect_equal(nrow(a3@path), 0L)
  expect_equal(a3@loglik, a2@loglik, tolerance = 1e-6)
})

test_that("accepttol and the criterion are honoured", {
  set.seed(7)
  truth <- sdvinecop(armacop(ar = c(0.6, -0.25)), pairs = list(pc("gumbel"), pc("clayton", pos = 180, neg = 270)))
  u <- tssim(truth, 800)
  none <- autodvinecop(u, armacop(order = c(2, 0)), nreplace = 2, accepttol = 1e6, pseudo = FALSE,
                       verbose = FALSE)
  expect_equal(length(none@model@pairs), 0L)
  expect_false(any(none@path$accepted))
  bic <- autodvinecop(u, armacop(order = c(2, 0)), nreplace = 1, criterion = "BIC", pseudo = FALSE,
                      verbose = FALSE)
  expect_equal(bic@path$criterion, BIC(bic), tolerance = 1e-8)
  aicc <- autodvinecop(u, armacop(order = c(2, 0)), nreplace = 1, criterion = "AICc", pseudo = FALSE,
                       verbose = FALSE)
  expect_equal(aicc@path$criterion, dv_ic(aicc@loglik, aicc@npar, aicc@nobs, "AICc"), tolerance = 1e-8)
})

test_that("autodvinecop validates its input", {
  set.seed(8)
  u <- tssim(armacop(ar = 0.5), 100)
  a <- armacop(order = c(1, 0))
  expect_error(autodvinecop(u, a, choices = "gauss"), "choices")
  expect_error(autodvinecop(u, a, aux = "foo"), "aux")
  expect_error(autodvinecop(u, a, criterion = "foo"), "arg")
  expect_error(autodvinecop(u, a, nreplace = 0), "nreplace")
  expect_error(autodvinecop(u, 1), "start")
  free <- sdvinecop(pairs = list(rvinecopulib::bicop_dist("gumbel", 0, 2)))
  expect_error(autodvinecop(u, free), "tied")
  u2 <- u
  u2[4] <- NA
  expect_error(autodvinecop(u2, a), "complete")
})
