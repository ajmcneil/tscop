test_that("pc specifications are validated", {
  expect_s4_class(pc("gumbel"), "pc")
  expect_equal(pc("Gaussian")@family, "gauss")
  expect_equal(pc("clayton", pos = 180, neg = 270)@pos, 180)
  expect_equal(c(pc("frank", pos = 180)@pos, pc("t", aux = 5)@neg), c(0, 0))
  expect_true(is.na(pc("gumbel")@aux))
  expect_equal(pc("t", aux = 6)@aux, 6)
  expect_error(pc("foo"), "family")
  expect_error(pc("gumbel", pos = 90), "pos")
  expect_error(pc("gumbel", neg = 180), "neg")
  expect_error(pc("gumbel", aux = 2), "auxiliary")
  expect_error(pc("t", aux = 1.5), "degrees")
  expect_error(pc("bb1", aux = 0.5), "BB1")
  expect_output(show(pc("gumbel")), "gumbel")
})

test_that("tied models are constructed and the effective lag is stored", {
  m <- sdvinecop(armacop(ar = 0.5, ma = 0.4), pairs = list(pc("gumbel"), pc("clayton")))
  expect_s4_class(m, "sdvinecop")
  expect_equal(m@mode, "tied")
  expect_equal(m@tail@family, "gauss")
  e <- arma_coef_expanded(m@kpacf)
  expect_equal(m@K, dv_effective_lag(e$ar, e$ma, 2L)$K)
  expect_equal(m@Klimit, "tol")
  # pure ARMA copula
  m0 <- sdvinecop(armacop(ar = 0.5, ma = 0.4))
  expect_equal(length(m0@pairs), 0)
  # pure AR: K is the AR order
  m1 <- sdvinecop(armacop(ar = c(0.5, 0.2)), pairs = list(pc("t", aux = 6)))
  expect_equal(m1@K, 2L)
  expect_equal(m1@Klimit, "ar")
  # floor at number of explicit pairs
  m2 <- sdvinecop(armacop(ar = 0.5), pairs = list(pc("gumbel"), pc("clayton"), pc("joe")))
  expect_equal(m2@K, 3L)
  expect_equal(m2@Klimit, "floor")
  # maxlag
  m3 <- sdvinecop(armacop(ar = 0.9, ma = 0.5), maxlag = 5)
  expect_equal(c(m3@K, m3@Klimit), c("5", "maxlag"), ignore_attr = TRUE)
})

test_that("tails are restricted to families that tend to independence", {
  for (f in c("gumbel", "clayton", "joe", "frank", "gauss")) {
    expect_s4_class(sdvinecop(armacop(ar = 0.5, ma = 0.4), tail = f), "sdvinecop")
  }
  expect_error(sdvinecop(armacop(ar = 0.5, ma = 0.4), tail = "t"), "tail family")
  expect_error(sdvinecop(armacop(ar = 0.5, ma = 0.4), tail = pc("bb1", aux = 1.5)), "tail family")
  m <- sdvinecop(armacop(ar = 0.5, ma = 0.4), tail = pc("clayton", pos = 180, neg = 270))
  expect_equal(m@tail@pos, 180)
})

test_that("construction errors are informative", {
  expect_error(sdvinecop(), "needs 'kpacf'")
  expect_error(sdvinecop(pairs = list(pc("gumbel"))), "needs 'kpacf'")
  expect_error(sdvinecop(kpacf = 1), "Gaussian copula process")
  expect_error(sdvinecop(armacop(ar = 0.5), pairs = list(pc("gumbel"), rvinecopulib::bicop_dist("gauss", 0, 0.3))), "mix")
  expect_error(sdvinecop(armacop(ar = 0.5), pairs = list(1)), "pc\\(\\)")
  expect_error(sdvinecop(armacop(ar = 0.5), tautol = 0), "tautol")
  expect_error(sdvinecop(armacop(ar = 0.5), maxlag = 0), "maxlag")
})

test_that("a model structure has an undetermined effective lag", {
  m <- sdvinecop(armacop(order = c(1, 1)), pairs = list(pc("gumbel")))
  expect_true(is.na(m@K))
  expect_output(show(m), "not yet determined")
  expect_error(dv_resolve(m), "unspecified")
})

test_that("free-mode models hold bicop_dist objects", {
  b <- list(rvinecopulib::bicop_dist("gumbel", 90, 2), rvinecopulib::bicop_dist("clayton", 0, 1))
  m <- sdvinecop(pairs = b)
  expect_equal(m@mode, "free")
  expect_equal(m@K, 2L)
  expect_null(m@tail)
  pr <- dv_resolve(m)
  expect_equal(pr[[1]]$rotation, 90)
  # Gaussian ARMA tail beyond the explicit copulas
  m2 <- sdvinecop(armacop(ar = 0.5, ma = 0.4), pairs = b)
  expect_equal(m2@mode, "free")
  expect_gt(m2@K, 2L)
  pr2 <- dv_resolve(m2)
  alpha <- arma_pacf(0.5, 0.4, m2@K)
  expect_equal(pr2[[5]]$parameters, alpha[5])
  expect_error(sdvinecop(armacop(ar = 0.5), pairs = b, tail = "gumbel"), "Gaussian tail")
  expect_output(show(m), "free mode")
})

test_that("tied resolution calibrates every lag to the ARMA kpacf", {
  m <- sdvinecop(armacop(ar = 0.5, ma = 0.4),
                pairs = list(pc("gumbel"), pc("clayton", pos = 180, neg = 270), pc("t", aux = 6)),
                tail = "frank")
  pr <- dv_resolve(m)
  expect_length(pr, m@K)
  tau <- dv_kpacf_tau(m@kpacf, m@K)
  for (j in seq_len(m@K)) expect_lt(abs(dv_pair_ktau(pr[[j]]) - tau[j]), 1e-7)
  expect_equal(pr[[3]]$parameters[2], 6)
  expect_equal(pr[[5]]$family, "frank")
})

test_that("negative lags use the negative rotations", {
  m <- sdvinecop(armacop(ar = -0.5), pairs = list(pc("gumbel", neg = 270)))
  expect_equal(dv_resolve(m)[[1]]$rotation, 270)
  expect_output(show(m), "gumbel270")
  m2 <- sdvinecop(armacop(ar = 0.5), pairs = list(pc("gumbel", pos = 180)))
  expect_output(show(m2), "gumbel180")
})

test_that("sdvinecop without substitutions reproduces armacop", {
  set.seed(5)
  for (cs in list(list(ar = 0.5, ma = 0.4), list(ar = c(0.5, 0.3)), list(ma = 0.6))) {
    arma <- do.call(armacop, cs)
    e <- arma_coef_expanded(arma)
    u <- pnorm(sim_arma_scores(40, e$ar, e$ma))
    m <- sdvinecop(arma, tautol = 1e-14)
    expect_equal(dv_model_loglik(m, u), arma_copula_loglik(u, e$ar, e$ma), tolerance = 1e-8)
    m2 <- sdvinecop(arma, pairs = lapply(seq_len(2), function(i) pc("gauss")), tautol = 1e-14)
    expect_equal(dv_model_loglik(m2, u), arma_copula_loglik(u, e$ar, e$ma), tolerance = 1e-8)
  }
})

test_that("truncation error is small at the default tolerance", {
  set.seed(6)
  arma <- armacop(ar = 0.7, ma = 0.4)
  e <- arma_coef_expanded(arma)
  u <- pnorm(sim_arma_scores(300, e$ar, e$ma))
  exact <- arma_copula_loglik(u, e$ar, e$ma)
  m <- sdvinecop(arma)
  expect_lt(abs(dv_model_loglik(m, u) - exact), 1e-2)
})

test_that("infeasible models give -Inf log-likelihood", {
  m <- sdvinecop(armacop(ar = 0.3), pairs = list(pc("bb1", aux = 3)))
  expect_equal(dv_model_loglik(m, runif(50)), -Inf)
})

test_that("coef lists ARMA coefficients then auxiliary parameters", {
  m <- sdvinecop(armacop(ar = 0.5, ma = 0.4),
                pairs = list(pc("gumbel"), pc("t", aux = 6), pc("bb1", aux = 1.4)))
  expect_equal(names(coef(m)), c("ar1", "ma1", "nu2", "delta3"))
  expect_equal(unname(coef(m)), c(0.5, 0.4, 6, 1.4))
})

test_that("show describes the model", {
  m <- sdvinecop(armacop(ar = 0.5, ma = 0.4), pairs = list(pc("gumbel"), pc("t", aux = 6)))
  expect_output(show(m), "tied mode")
  expect_output(show(m), "gumbel0, t")
  expect_output(show(m), "effective lag K")
})

test_that("log-likelihood agrees with the legacy sdvinecopula objective", {
  skip_if_no_oracle()
  set.seed(7)
  u <- pnorm(sim_arma_scores(200, 0.5, 0.3))
  check <- function(ar, ma, fam, posrot, negrot, auxpar, pcs) {
    mod <- sdvinecop(armacop(ar = ar, ma = ma), pairs = pcs)
    old <- tscopula::sdvinecopula(pars = list(ar = ar, ma = ma), family = fam, posrot = posrot,
                                  negrot = negrot, auxpar = auxpar, tautol = 1e-4)
    th <- unlist(old@pars)
    names(th) <- sub("^[^.]*\\.", "", names(th))
    expect_equal(dv_model_loglik(mod, u),
                 -tscopula:::sdvinecopula_objective(th, old@modelspec, u), tolerance = 1e-8)
  }
  check(0.5, 0.3, c("gumbel", "clayton", "t"), c(0, 180, 0), c(90, 270, 0), 6,
        list(pc("gumbel"), pc("clayton", pos = 180, neg = 270), pc("t", aux = 6)))
  check(-0.4, 0.3, c("gumbel", "clayton", "joe"), c(0, 180, 0), c(90, 270, 90), NA,
        list(pc("gumbel"), pc("clayton", pos = 180, neg = 270), pc("joe")))
  check(0.6, -0.5, c("frank", "gumbel"), c(0, 0), c(0, 90), NA, list(pc("frank"), pc("gumbel")))
})

test_that("rep repeats a pc specification and pairs accepts a single pc", {
  r <- rep(pc("ivt", base = "ast"), 3)
  expect_type(r, "list")
  expect_length(r, 3L)
  expect_true(all(vapply(r, is, logical(1), "pc")))
  expect_equal(length(rep(pc("gumbel"), times = 2, each = 2)), 4L)
  expect_equal(length(rep(pc("gumbel"), length.out = 5)), 5L)
  kp <- armacop(ar = 0.5, ma = 0.4)
  expect_equal(coef(sdvinecop(kp, pairs = rep(pc("gumbel"), 2))),
               coef(sdvinecop(kp, pairs = list(pc("gumbel"), pc("gumbel")))))
  m1 <- sdvinecop(kp, pairs = pc("gumbel"))
  expect_equal(length(m1@pairs), 1L)
  mixed <- sdvinecop(kp, pairs = c(rep(pc("gumbel"), 2), rep(pc("gauss"), 2)))
  expect_equal(length(mixed@pairs), 4L)
  expect_error(sdvinecop(kp, pairs = "gumbel"), "list")
})
