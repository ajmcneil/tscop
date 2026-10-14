# The compiled simulation and Rosenblatt kernels against the R implementations.

sim_models <- function() {
  kp <- arfactcop(lambda = 0.95, w = 0.5)
  arm <- armacop(ar = c(0.5, -0.3))
  list(
    plain_joe = sdvinecop(kp, tail = pc("joe"), maxlag = 12, tautol = 1e-6),
    plain_gumbel = sdvinecop(kp, tail = pc("gumbel"), maxlag = 12, tautol = 1e-6),
    plain_clayton = sdvinecop(kp, tail = pc("clayton"), maxlag = 12, tautol = 1e-6),
    plain_frank = sdvinecop(kp, tail = pc("frank"), maxlag = 12, tautol = 1e-6),
    plain_gauss = sdvinecop(kp, tail = pc("gauss"), maxlag = 12, tautol = 1e-6),
    rot_0_90 = sdvinecop(arm, pairs = list(pc("gumbel"), pc("clayton"), pc("joe"), pc("frank")), tail = "gauss", maxlag = 8),
    rot_180_270 = sdvinecop(arm, pairs = list(pc("gumbel", pos = 180, neg = 270), pc("clayton", pos = 180, neg = 270),
                                             pc("joe", pos = 180, neg = 270)), tail = pc("gumbel", pos = 180, neg = 270), maxlag = 8),
    ivt_joe = sdvinecop(kp, tail = pc("ivt", base = "joe"), vt = c(0.4, 0.6), maxlag = 10, tautol = 1e-6),
    ivt_gumbel = sdvinecop(kp, tail = pc("ivt", base = "gumbel"), vt = c(0.4, 0.6), maxlag = 10, tautol = 1e-6),
    ivt_clayton180 = sdvinecop(kp, tail = pc("ivt", base = "clayton180"), vt = c(0.4, 0.6), maxlag = 10, tautol = 1e-6),
    ivt_frank = sdvinecop(kp, tail = pc("ivt", base = "frank"), vt = c(0.4, 0.6), maxlag = 10, tautol = 1e-6),
    ivt_gauss = sdvinecop(kp, tail = pc("ivt", base = "gauss"), vt = c(0.4, 0.6), maxlag = 10, tautol = 1e-6),
    mixed = sdvinecop(kp, pairs = c(rep(pc("ivt", base = "joe"), 3), list(pc("gumbel")), rep(pc("ivt", base = "frank"), 2)),
                      tail = pc("ivt", base = "gauss"), vt = c(0.4, 0.6), maxlag = 10)
  )
}

test_that("the C++ Rosenblatt transform equals the R one for every supported pair copula", {
  set.seed(1)
  for (nm in names(sim_models())) {
    m <- sim_models()[[nm]]
    m@K <- as.integer(min(m@maxlag, 12))
    pr <- dv_resolve(m)
    u <- runif(300)
    expect_equal(dv_rosenblatt(pr, u), dv_rosenblatt_r(pr, u), tolerance = 1e-8, info = nm)
  }
})

test_that("the C++ simulation inverts the Rosenblatt transform and agrees with the R simulation", {
  set.seed(2)
  for (nm in names(sim_models())) {
    m <- sim_models()[[nm]]
    m@K <- as.integer(min(m@maxlag, 12))
    pr <- dv_resolve(m)
    expect_false(is.null(dv_encode(pr)))
    w <- matrix(runif(250 * 2), 250, 2)
    u <- dv_simulate(pr, w)
    expect_true(all(u > 0 & u < 1), info = nm)
    # the compiled inverse is exact: the Rosenblatt transform returns the innovations
    for (j in 1:2) expect_equal(dv_rosenblatt(pr, u[, j]), w[, j], tolerance = 1e-9, info = nm)
    # rvinecopulib's numerical inverses are accurate to only about 1e-5 (Joe)
    expect_equal(u, dv_simulate_r(pr, w), tolerance = 1e-3, info = nm)
  }
})

test_that("simulation and the Rosenblatt transform fall back to R for unsupported families", {
  set.seed(3)
  m <- sdvinecop(arfactcop(lambda = 0.9, w = 0.5), pairs = list(pc("t", aux = 6), pc("gumbel")), tail = "gauss", maxlag = 5)
  m@K <- 5L
  pr <- dv_resolve(m)
  expect_null(dv_encode(pr))
  w <- matrix(runif(100), 100, 1)
  expect_equal(dv_simulate(pr, w), dv_simulate_r(pr, w))
  u <- as.numeric(dv_simulate(pr, w))
  expect_equal(dv_rosenblatt(pr, u), dv_rosenblatt_r(pr, u))
  old <- options(tscop.cpp = FALSE)
  on.exit(options(old))
  mj <- sdvinecop(arfactcop(lambda = 0.9, w = 0.5), tail = pc("joe"), maxlag = 5)
  mj@K <- 5L
  expect_equal(dv_simulate(dv_resolve(mj), w), dv_simulate_r(dv_resolve(mj), w))
})

test_that("tssim, resid and the bootstrap use the kernels and give consistent results", {
  set.seed(4)
  m <- sdvinecop(arfactcop(lambda = 0.95, w = 0.5), tail = pc("ivt", base = "joe"), vt = c(0.4, 0.6), maxlag = 15, tautol = 1e-4)
  u <- tssim(m, 500)
  expect_true(all(u > 0 & u < 1))
  r <- resid(m, u)
  expect_gt(stats::ks.test(r, "punif")$p.value, 0.001)
  expect_lt(abs(stats::cor(r[-1], r[-length(r)])), 0.15)
})

# ---- predictive distribution -------------------------------------------------

test_that("the C++ predictive distribution equals the R one for every supported pair copula", {
  set.seed(5)
  x <- c(0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99)
  for (nm in names(sim_models())) {
    m <- sim_models()[[nm]]
    m@K <- as.integer(min(m@maxlag, 12))
    pr <- dv_resolve(m)
    u <- tssim(m, 150)
    for (type in c("df", "qf", "dens")) {
      # the rotated model has Frank pairs with tau near zero, where rvinecopulib's density is inaccurate
      tol <- if (nm == "rot_0_90" && type == "dens") 1e-3 else 1e-8
      expect_equal(dv_predictive(pr, u, x, type), dv_predictive_r(pr, u, x, type), tolerance = tol, info = paste(nm, type))
      # fewer past values than the order of the model
      expect_equal(dv_predictive(pr, u[1:3], x, type), dv_predictive_r(pr, u[1:3], x, type), tolerance = tol, info = paste(nm, type, "short"))
    }
    expect_equal(dv_predictive(pr, u, dv_predictive(pr, u, x, "df"), "qf"), x, tolerance = 1e-8, info = nm)
  }
})

test_that("the predictive density integrates to one and matches the distribution function", {
  set.seed(6)
  m <- sdvinecop(arfactcop(lambda = 0.95, w = 0.5), tail = pc("ivt", base = "joe"), vt = c(0.4, 0.6), maxlag = 10, tautol = 1e-6)
  m@K <- 10L
  pr <- dv_resolve(m)
  u <- tssim(m, 100)
  g <- seq(0, 1, length.out = 4001)[-c(1, 4001)]
  d <- dv_predictive(pr, u, g, "dens")
  expect_equal(sum(d) / length(g), 1, tolerance = 5e-3)
  Fx <- dv_predictive(pr, u, g, "df")
  expect_equal(diff(Fx) / diff(g), 0.5 * (d[-1] + d[-length(d)]), tolerance = 5e-2)
})

test_that("predict for d-vine and v-transform models uses the kernel and agrees with R", {
  set.seed(7)
  m <- sdvinecop(arfactcop(lambda = 0.95, w = 0.5), tail = pc("joe"), maxlag = 8, tautol = 1e-6)
  u <- tssim(m, 120)
  p1 <- predict(m, u, type = "qf", x = c(0.05, 0.5, 0.95))
  old <- options(tscop.cpp = FALSE)
  on.exit(options(old))
  p2 <- predict(m, u, type = "qf", x = c(0.05, 0.5, 0.95))
  expect_equal(p1, p2, tolerance = 1e-8)
})

test_that("simulation, Rosenblatt and predictive kernels cover BB1 pairs with rotations", {
  set.seed(26)
  for (nm in names(bb1_models())) {
    m <- bb1_models()[[nm]]
    m@K <- 10L
    pr <- dv_resolve(m)
    expect_false(is.null(dv_encode(pr)))
    w <- matrix(runif(200 * 2), 200, 2)
    u <- dv_simulate(pr, w)
    for (j in 1:2) expect_equal(dv_rosenblatt(pr, u[, j]), w[, j], tolerance = 1e-9, info = nm)
    expect_equal(u, dv_simulate_r(pr, w), tolerance = 1e-3, info = nm)
    x <- c(0.1, 0.5, 0.9)
    for (type in c("df", "qf", "dens")) {
      expect_equal(dv_predictive(pr, u[, 1], x, type), dv_predictive_r(pr, u[, 1], x, type), tolerance = 1e-7, info = paste(nm, type))
    }
  }
})

test_that("simulation, Rosenblatt and predictive kernels cover ast pairs, plain and as an IVT base", {
  set.seed(34)
  for (nm in names(ast_models())) {
    m <- ast_models()[[nm]]
    m@K <- 10L
    pr <- dv_resolve(m)
    expect_false(is.null(dv_encode(pr, ast = TRUE)))
    w <- matrix(runif(200 * 2), 200, 2)
    u <- dv_simulate(pr, w)
    for (j in 1:2) expect_equal(dv_rosenblatt(pr, u[, j]), w[, j], tolerance = 1e-9, info = nm)
    expect_equal(u, dv_simulate_r(pr, w), tolerance = 1e-3, info = nm)
    x <- c(0.05, 0.5, 0.95)
    for (type in c("df", "qf", "dens")) {
      expect_equal(dv_predictive(pr, u[, 1], x, type), dv_predictive_r(pr, u[, 1], x, type), tolerance = 1e-7, info = paste(nm, type))
    }
  }
})
