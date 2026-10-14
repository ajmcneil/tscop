test_that("the O(m log m) Kendall's tau equals stats::cor, with and without ties", {
  set.seed(1)
  for (n in c(2, 5, 50, 700)) {
    x <- rnorm(n)
    y <- 0.5 * x + rnorm(n)
    expect_equal(kendall_tau_cpp(x, y), stats::cor(x, y, method = "kendall"), tolerance = 1e-12)
  }
  x <- sample(1:6, 300, TRUE)
  y <- sample(1:5, 300, TRUE) + x %/% 2
  expect_equal(kendall_tau_cpp(x, y), stats::cor(x, y, method = "kendall"), tolerance = 1e-12)
  expect_equal(c(kendall_tau_cpp(1:10, 1:10), kendall_tau_cpp(1:10, 10:1)), c(1, -1))
  expect_true(is.nan(kendall_tau_cpp(1:10, rep(1, 10))))
  expect_error(kendall_tau_cpp(1:3, 1:4), "same length")
})

test_that("the compiled level taus equal the R reference for every pair copula family", {
  set.seed(2)
  models <- c(ast_models(), t_models(), bb1_models())
  for (nm in names(models)) {
    m <- models[[nm]]
    m@K <- 10L
    pr <- dv_resolve(m)
    u <- tssim(m, 600)
    expect_equal(dv_level_tau(pr, u), dv_level_tau_r(pr, u)$tau, tolerance = 1e-3, info = nm)
  }
  m <- sdvinecop(arfactcop(lambda = 0.9, w = 0.5), tail = pc("ivt", base = "gumbel"), vt = c(0.3, 0.7), maxlag = 10, tautol = 1e-6)
  m@K <- 10L
  pr <- dv_resolve(m)
  u <- tssim(m, 800)
  expect_equal(dv_level_tau(pr, u), dv_level_tau_r(pr, u)$tau, tolerance = 1e-3)
})

test_that("ktauhat reproduces the generalised-lag algorithm computed with rvinecopulib (independent oracle)", {
  set.seed(3)
  m <- armacop(ar = c(0.5, 0.2), ma = 0.3)
  u <- tssim(m, 800)
  r <- ktauhat(m, u, lagmax = 6)
  pacf <- gc_pacf(m, 6)
  d <- cbind(u[-800], u[-1])
  ref <- numeric(6)
  ref[1] <- stats::cor(d, method = "kendall")[1, 2]
  for (i in 1:5) {
    n <- nrow(d)
    mod <- rvinecopulib::bicop_dist("gaussian", 0, pacf[i])
    d <- cbind(rvinecopulib::hbicop(d[1:(n - 1), ], mod, cond_var = 2), rvinecopulib::hbicop(d[2:n, ], mod, cond_var = 1))
    ref[i + 1] <- stats::cor(d, method = "kendall")[1, 2]
  }
  expect_equal(r$tau_data, ref, tolerance = 1e-6)
  expect_equal(r$tau_model, ktau(m, 6))
  expect_equal(r$n, 800 - 1:6)
  expect_equal(names(r), c("lag", "tau_data", "tau_model", "n"))
})

test_that("on data from the model the empirical Kendall partial autocorrelations match the model's kpacf", {
  skip_on_cran()
  set.seed(4)
  kp <- arfactcop(lambda = 0.9, w = 0.5)
  models <- list(
    armacop = armacop(ar = 0.6, ma = 0.3),
    arfactcop = kp,
    dvine = sdvinecop(kp, pairs = list(pc("gumbel"), pc("t", aux = 6)), tail = pc("joe"), maxlag = 12, tautol = 1e-3),
    ivt_ast = sdvinecop(kp, tail = pc("ivt", base = "ast"), vt = c(0.45, 0.55), maxlag = 12, tautol = 1e-3),
    vtscop = vtscop(kp, udp::v2p(0.45, 0.8))
  )
  for (nm in names(models)) {
    m <- models[[nm]]
    u <- tssim(m, 20000)
    r <- ktauhat(m, u, lagmax = 6)
    expect_lt(max(abs(r$tau_data - r$tau_model)), 0.03, label = nm)
  }
})

test_that("the simulation band covers the empirical values of model-true data and misses a wrong model", {
  skip_on_cran()
  set.seed(5)
  kp <- arfactcop(lambda = 0.9, w = 0.5)
  m <- sdvinecop(kp, pairs = list(pc("gumbel")), maxlag = 12, tautol = 1e-3)
  u <- tssim(m, 3000)
  r <- ktauhat(m, u, lagmax = 8, B = 80)
  expect_true(all(c("lo", "hi") %in% names(r)))
  expect_true(all(r$lo < r$hi))
  expect_gte(mean(r$tau_data >= r$lo & r$tau_data <= r$hi), 0.75)
  expect_equal(attr(r, "level"), 0.95)
  wrong <- sdvinecop(arfactcop(lambda = 0.5, w = 0.5), pairs = list(pc("gumbel")), maxlag = 12, tautol = 1e-3)
  rw <- ktauhat(wrong, u, lagmax = 8, B = 80)
  expect_lt(mean(rw$tau_data >= rw$lo & rw$tau_data <= rw$hi), 0.75)
})

test_that("the pairs argument returns the first level samples and is capped at 9", {
  set.seed(6)
  m <- sdvinecop(arfactcop(lambda = 0.9, w = 0.5), pairs = list(pc("gumbel")), maxlag = 12, tautol = 1e-3)
  u <- tssim(m, 500)
  expect_null(attr(ktauhat(m, u, lagmax = 6), "pairs"))
  expect_null(attr(ktauhat(m, u, lagmax = 6, pairs = NULL), "pairs"))
  r <- ktauhat(m, u, lagmax = 6, pairs = 3)
  p <- attr(r, "pairs")
  expect_length(p, 3L)
  expect_equal(dim(p[[1]]), c(499L, 2L))
  expect_equal(dim(p[[3]]), c(497L, 2L))
  expect_equal(unname(p[[1]]), cbind(u[-500], u[-1]), ignore_attr = TRUE)
  expect_true(all(p[[3]] > 0 & p[[3]] < 1))
  # the tau of the returned pairs is the reported tau (for a model without v-transforms)
  expect_equal(stats::cor(p[[3]], method = "kendall")[1, 2], r$tau_data[3], tolerance = 1e-9)
  expect_warning(r9 <- ktauhat(m, u, lagmax = 12, pairs = 15), "capped at 9")
  expect_length(attr(r9, "pairs"), 9L)
  expect_length(attr(ktauhat(m, u, lagmax = 2, pairs = 9), "pairs"), 2L)
})

test_that("ktauhat handles fitted objects, full models, a changed backend and bad input", {
  set.seed(7)
  kp <- arfactcop(lambda = 0.9, w = 0.5)
  tm <- tsmod(marg("st", c(mu = 1, sigma = 2, df = 6)), sdvinecop(kp, pairs = list(pc("gumbel")), maxlag = 12, tautol = 1e-3))
  x <- tssim(tm, 800)
  f <- tsfit(x, tsmod("st", sdvinecop(arfactcop(order = 1), pairs = list(pc("gumbel")), maxlag = 12, tautol = 1e-3)))
  r_fit <- ktauhat(f, lagmax = 5)
  r_mod <- ktauhat(f@model, x, lagmax = 5)
  expect_equal(r_fit, r_mod)
  expect_equal(names(r_fit), c("lag", "tau_data", "tau_model", "n"))
  u <- tssim(f@model@copula, 800)
  cf <- tsfit(u, sdvinecop(arfactcop(order = 1), pairs = list(pc("gumbel")), maxlag = 12, tautol = 1e-3), pseudo = FALSE)
  expect_equal(ktauhat(cf, lagmax = 5), ktauhat(cf@model, cf@data, lagmax = 5))
  # the R backend gives the same table
  old <- options(tscop.cpp = FALSE)
  on.exit(options(old))
  expect_equal(ktauhat(f, lagmax = 5)$tau_data, r_fit$tau_data, tolerance = 1e-4)
  options(old)
  # input checks
  expect_error(ktauhat(f@model@copula), "'data' is required")
  u2 <- u
  u2[3] <- NA
  expect_error(ktauhat(f@model@copula, u2), "without missing values")
  expect_error(ktauhat(f@model@copula, u, lagmax = 0), "positive whole number")
  expect_equal(nrow(ktauhat(f@model@copula, u[1:40], lagmax = 100)), 28L)       # capped at n - 12
  expect_error(ktauhat(f@model@copula, u[1:12]), "too few")
  expect_error(ktauhat(armacop(order = c(1, 1)), u), "unspecified")
})

test_that("the simulation band does not depend on the number of cores and leaves the random number stream as the caller's draw", {
  skip_on_cran()
  skip_on_os("windows")
  kp <- arfactcop(lambda = 0.9, w = 0.5)
  m <- sdvinecop(kp, pairs = list(pc("gumbel")), maxlag = 12, tautol = 1e-3)
  set.seed(8)
  u <- tssim(m, 600)
  set.seed(9)
  r1 <- ktauhat(m, u, lagmax = 5, B = 12)
  set.seed(9)
  r2 <- ktauhat(m, u, lagmax = 5, B = 12, cores = 2)
  expect_equal(r1, r2)
  # the caller's stream advances by one draw of B seeds, whatever the number of cores
  set.seed(9)
  invisible(ktauhat(m, u, lagmax = 5, B = 12))
  after <- runif(1)
  set.seed(9)
  invisible(sample.int(.Machine$integer.max, 12))
  expect_equal(after, runif(1))
  # the other classes
  set.seed(10)
  expect_true(all(c("lo", "hi") %in% names(ktauhat(armacop(ar = 0.6), tssim(armacop(ar = 0.6), 400), lagmax = 4, B = 10))))
  v <- vtscop(kp, udp::v2p(0.45, 0.8))
  expect_true(all(c("lo", "hi") %in% names(ktauhat(v, tssim(v, 600), lagmax = 4, B = 10))))
  tm <- tsmod(marg("st", c(mu = 0, sigma = 1, df = 5)), m)
  expect_true(all(c("lo", "hi") %in% names(ktauhat(tm, tssim(tm, 600), lagmax = 4, B = 10))))
})
