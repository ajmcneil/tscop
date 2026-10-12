# ---- the generic rotation layer matches rvinecopulib's rotated families ----

rv_closures <- function(family, parameters) {
  list(
    d = function(a, b) as.numeric(rvinecopulib::dbicop(cbind(a, b), family, 0, parameters)),
    h = function(a, b, cv) as.numeric(rvinecopulib::hbicop(cbind(a, b), cv, family, 0, parameters)),
    hinv = function(a, p) as.numeric(rvinecopulib::hbicop(cbind(a, p), 1L, family, 0, parameters, inverse = TRUE))
  )
}

test_that("dv_rotate reproduces rvinecopulib rotations (density, both h-functions, inverse)", {
  set.seed(1)
  a <- runif(30, 0.02, 0.98)
  b <- runif(30, 0.02, 0.98)
  for (cs in list(list("gumbel", 1.8), list("clayton", 2.2), list("joe", 2.5), list("frank", 4))) {
    base <- rv_closures(cs[[1]], cs[[2]])
    for (rot in c(0, 90, 180, 270)) {
      if (cs[[1]] == "frank" && rot != 0) next
      wrapped <- dv_rotate(base, rot)
      direct <- list(
        d = as.numeric(rvinecopulib::dbicop(cbind(a, b), cs[[1]], rot, cs[[2]])),
        h1 = as.numeric(rvinecopulib::hbicop(cbind(a, b), 1, cs[[1]], rot, cs[[2]])),
        h2 = as.numeric(rvinecopulib::hbicop(cbind(a, b), 2, cs[[1]], rot, cs[[2]])),
        hinv = as.numeric(rvinecopulib::hbicop(cbind(a, b), 1, cs[[1]], rot, cs[[2]], inverse = TRUE))
      )
      expect_equal(wrapped$d(a, b), direct$d, tolerance = 1e-8, info = paste(cs[[1]], rot))
      expect_equal(wrapped$h(a, b, 1L), direct$h1, tolerance = 1e-8, info = paste(cs[[1]], rot))
      expect_equal(wrapped$h(a, b, 2L), direct$h2, tolerance = 1e-8, info = paste(cs[[1]], rot))
      expect_equal(wrapped$hinv(a, b), direct$hinv, tolerance = 1e-7, info = paste(cs[[1]], rot))
    }
  }
  expect_error(dv_rotate(rv_closures("gumbel", 2), 45), "rotation")
})

# ---- the ast pair copula -------------------------------------------------

test_that("ast pair functions are mutually consistent for every rotation", {
  set.seed(2)
  a <- runif(15, 0.1, 0.9)
  b <- runif(15, 0.1, 0.9)
  for (rot in c(0, 90, 180, 270)) {
    pr <- list(family = "ast", rotation = rot, parameters = 1.5)
    e <- 1e-6
    # density is the b-derivative of h(b | a) and the a-derivative of h(a | b)
    d1 <- (dv_h(pr, a, b + e, 1L) - dv_h(pr, a, b - e, 1L)) / (2 * e)
    d2 <- (dv_h(pr, a + e, b, 2L) - dv_h(pr, a - e, b, 2L)) / (2 * e)
    expect_equal(d1, dv_dens(pr, a, b), tolerance = 1e-4, info = rot)
    expect_equal(d2, dv_dens(pr, a, b), tolerance = 1e-4, info = rot)
    # inverse
    expect_equal(dv_hinv(pr, a, dv_h(pr, a, b, 1L)), b, tolerance = 1e-6, info = rot)
    # h-functions are conditional distribution functions
    expect_true(all(dv_h(pr, a, b, 1L) > 0 & dv_h(pr, a, b, 1L) < 1))
  }
  # rotations 90 and 270 are negatively dependent, 0 and 180 positively
  expect_gt(dv_pair_ktau(list(family = "ast", rotation = 0, parameters = 1.5)), 0)
  expect_gt(dv_pair_ktau(list(family = "ast", rotation = 180, parameters = 1.5)), 0)
  expect_lt(dv_pair_ktau(list(family = "ast", rotation = 90, parameters = 1.5)), 0)
  expect_lt(dv_pair_ktau(list(family = "ast", rotation = 270, parameters = 1.5)), 0)
})

test_that("ast pairs are calibrated to Kendall's tau and rotated by sign", {
  for (tau in c(-0.8, -0.3, -0.02, 0.02, 0.3, 0.8)) {
    pr <- dv_resolve_vals("ast", 180, 270, NA, tau)
    expect_equal(pr$family, "ast")
    expect_equal(pr$rotation, if (tau >= 0) 180 else 270)
    expect_equal(dv_pair_ktau(pr), tau, tolerance = 1e-6)
  }
  expect_equal(dv_resolve_vals("ast", 0, 90, NA, 0)$family, "indep")
  expect_equal(dv_resolve_vals("ast", 0, 90, NA, 1e-12)$family, "indep")
  tail <- dv_resolve_tail(pc("ast", pos = 180, neg = 270), c(0.5, -0.1, 0, 0.01))
  expect_equal(vapply(tail, function(p) p$family, ""), c("ast", "ast", "indep", "ast"))
  expect_equal(vapply(tail, function(p) p$rotation, 0), c(180, 270, 0, 180))
  expect_equal(dv_pair_ktau(tail[[1]]), 0.5, tolerance = 1e-6)
})

# a plain loop through the pair-level functions (no z-scale shortcuts)
dv_loglik_ref2 <- function(u, pairs) {
  n <- length(u)
  K <- min(length(pairs), n - 1L)
  A <- u[-n]
  B <- u[-1]
  ll <- 0
  for (j in seq_len(K)) {
    ll <- ll + sum(log(dv_dens(pairs[[j]], A, B)))
    if (j < K) {
      m <- length(A)
      Anew <- dv_h(pairs[[j]], A[-m], B[-m], 2L)
      Bnew <- dv_h(pairs[[j]], A[-1], B[-1], 1L)
      A <- Anew
      B <- Bnew
    }
  }
  ll
}

test_that("the ast engine agrees with a plain loop over the udp functions", {
  set.seed(3)
  u <- runif(80)
  nus <- c(1.2, 2.5, 6)
  pairs <- lapply(nus, function(n) list(family = "ast", rotation = 0, parameters = n))
  ref <- 0
  v <- cbind(u[-80], u[-1])
  for (j in seq_along(pairs)) {
    ac <- udp::astcopula(nus[j])
    ref <- ref + sum(log(udp::dastcopula(v[, 1], v[, 2], ac)))
    if (j < length(pairs)) {
      m <- nrow(v)
      v <- cbind(udp::hastcopula(v[-m, 1], v[-m, 2], ac, cond_var = 2L),
                 udp::hastcopula(v[-1, 1], v[-1, 2], ac, cond_var = 1L))
    }
  }
  expect_equal(dv_loglik(u, pairs), ref, tolerance = 1e-9)
  # mixed with other families
  mixed <- list(pair("gumbel", 0, 1.5), list(family = "ast", rotation = 0, parameters = 2),
                pair("gauss", 0, 0.2), list(family = "ast", rotation = 180, parameters = 4))
  expect_equal(dv_loglik(u, mixed), dv_loglik_ref2(u, mixed), tolerance = 1e-9)
})

test_that("simulation and Rosenblatt are inverse for models with ast pairs", {
  set.seed(4)
  pairs <- list(list(family = "ast", rotation = 0, parameters = 1.5),
                list(family = "ast", rotation = 270, parameters = 3),
                pair("gumbel", 0, 1.3),
                list(family = "ast", rotation = 180, parameters = 8))
  w <- matrix(runif(25 * 3, 0.02, 0.98), 25, 3)
  u <- dv_simulate(pairs, w)
  for (r in 1:3) expect_equal(dv_rosenblatt(pairs, u[, r]), w[, r], tolerance = 1e-6)
})

# ---- the pc and sdvinecop interface -----------------------------------------

test_that("ast is a valid pair copula family and tail", {
  p <- pc("ast", pos = 180, neg = 270)
  expect_equal(c(p@family, p@pos, p@neg), c("ast", "180", "270"))
  expect_true(is.na(p@aux))
  expect_error(pc("ast", aux = 2), "auxiliary")
  expect_s4_class(sdvinecop(armacop(ar = 0.5, ma = 0.4), tail = "ast"), "sdvinecop")
  m <- sdvinecop(armacop(ar = 0.5, ma = 0.4), pairs = list(pc("ast"), pc("gumbel")), tail = pc("ast"))
  expect_output(show(m), "tail family: ast")
  pr <- dv_resolve(m)
  expect_equal(pr[[1]]$family, "ast")
  tau <- dv_kpacf_tau(m@kpacf, m@K)
  for (j in seq_len(m@K)) expect_lt(abs(dv_pair_ktau(pr[[j]]) - tau[j]), 1e-6)
})

test_that("a d-vine with ast pair copulas simulates, scores and refits", {
  set.seed(5)
  truth <- sdvinecop(armacop(ar = 0.6), pairs = list(pc("ast")))
  u <- tssim(truth, 600)
  expect_true(all(u > 0 & u < 1))
  r <- resid(truth, u)
  expect_gt(ks.test(r, "punif")$p.value, 0.001)
  f <- tsfit(u, sdvinecop(armacop(order = c(1, 0)), pairs = list(pc("ast"))), pseudo = FALSE)
  expect_equal(unname(coef(f)), 0.6, tolerance = 0.12)
  fg <- tsfit(u, sdvinecop(armacop(order = c(1, 0))), pseudo = FALSE)
  expect_gt(f@loglik, fg@loglik)
  # predictive distribution functions are consistent
  x <- c(0.1, 0.5, 0.9)
  df <- predict(truth, u, type = "df", x = x)
  expect_equal(predict(truth, u, type = "qf", x = df), x, tolerance = 1e-5)
})

test_that("autodvinecop can select the ast copula", {
  set.seed(6)
  truth <- sdvinecop(armacop(ar = 0.7), pairs = list(pc("ast")))
  u <- tssim(truth, 700)
  gauss <- tsfit(u, armacop(order = c(1, 0)), pseudo = FALSE)
  a <- autodvinecop(u, gauss, nreplace = 1, choices = "ast", pseudo = FALSE, verbose = FALSE)
  expect_equal(a@path$choice, "ast0")
  expect_true(a@path$accepted)
  expect_lt(AIC(a), AIC(tsfit(u, sdvinecop(armacop(order = c(1, 0))), pseudo = FALSE)) - 10)
  expect_error(autodvinecop(u, gauss, choices = "t"), "choices")
})
