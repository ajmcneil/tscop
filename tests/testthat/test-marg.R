mpars <- list(
  gauss = c(mu = 0.3, sigma = 1.7),
  laplace = c(mu = -0.2, sigma = 0.8),
  slaplace = c(mu = 0.1, sigma = 1.2, gamma = 0.7),
  st = c(mu = 0.5, sigma = 0.9, df = 4.5),
  sst = c(mu = -0.3, sigma = 1.1, df = 6, gamma = 1.4),
  jsu = c(mu = 0.2, sigma = 1.5, delta = 2.2, gamma = -0.6),
  nig = c(mu = 0.2, sigma = 1.5, shape = 1.2, skew = -0.5),
  hyp = c(mu = -0.1, sigma = 0.9, shape = 0.8, skew = 0.6),
  gh = c(mu = 0.1, sigma = 1.2, lambda = 1.5, shape = 0.8, skew = 0.4)
)

test_that("marg validates its arguments", {
  m <- marg("sst", c(mu = 0, sigma = 1, df = 5, gamma = 1))
  expect_s4_class(m, "marg")
  expect_equal(names(m@pars), c("mu", "sigma", "df", "gamma"))
  expect_true(anyNA(marg("st")@pars))
  expect_error(marg("foo"), "unknown margin family")
  expect_error(marg("st", c(nu = 3)), "named with parameters")
  expect_error(marg("st", c(mu = 0, sigma = -1, df = 4)), "sigma")
  expect_error(marg("st", c(mu = 0, sigma = 1, df = 0.5)), "df")
  expect_error(marg("st", fix = "mu"), "fixed parameters need values")
  expect_error(marg("st", c(mu = 0), fix = "foo"), "fix")
  expect_equal(unname(coef(marg("st", c(mu = 0, sigma = 1, df = 5), fix = "mu"))), c(1, 5))
  expect_error(tssim(marg("st"), 10), "unspecified")
  expect_output(show(m), "skew Student t")
})

test_that("density, distribution and quantile functions are mutually consistent", {
  for (fam in names(mpars)) {
    m <- marg(fam, mpars[[fam]])
    th <- mpars[[fam]]
    x <- seq(th[["mu"]] - 8 * th[["sigma"]], th[["mu"]] + 8 * th[["sigma"]], length.out = 41)
    # density integrates to one
    expect_equal(stats::integrate(function(v) mdens(m, v), -Inf, Inf, rel.tol = 1e-8)$value, 1,
                 tolerance = 1e-6, info = fam)
    # derivative of the df is the density
    e <- 1e-5
    num <- (mdf(m, x + e) - mdf(m, x - e)) / (2 * e)
    expect_equal(num, mdens(m, x), tolerance = 1e-5, info = fam)
    expect_equal(mdens(m, x, log = TRUE), log(mdens(m, x)), info = fam)
    # tails add up and the quantile function inverts the df
    expect_equal(mdf(m, x) + mdf(m, x, lower.tail = FALSE), rep(1, length(x)), tolerance = 1e-12, info = fam)
    p <- c(1e-10, 1e-4, 0.01, 0.2, 0.5, 0.8, 0.99, 1 - 1e-4)
    expect_equal(mdf(m, mqf(m, p)), p, tolerance = 1e-9, info = fam)
    expect_equal(mdf(m, mqf(m, p, lower.tail = FALSE), lower.tail = FALSE), p, tolerance = 1e-9, info = fam)
    expect_equal(mqf(m, p[-1], lower.tail = FALSE), mqf(m, 1 - p[-1]), tolerance = 1e-8, info = fam)  # 1 - 1e-10 loses digits
    expect_true(all(diff(mdf(m, x)) > 0), info = fam)
  }
  # extreme upper tails keep their accuracy
  m <- marg("sst", mpars$sst)
  expect_gt(mdf(m, 1e3, lower.tail = FALSE), 0)
  expect_lt(mdf(m, 1e3, lower.tail = FALSE), 1e-6)
})

test_that("nested families coincide", {
  x <- seq(-4, 4, length.out = 21)
  expect_equal(mdens(marg("sst", c(mu = 0.2, sigma = 1.3, df = 5, gamma = 1)), x),
               mdens(marg("st", c(mu = 0.2, sigma = 1.3, df = 5)), x), tolerance = 1e-12)
  expect_equal(mdf(marg("slaplace", c(mu = 0.2, sigma = 1.3, gamma = 1)), x),
               mdf(marg("laplace", c(mu = 0.2, sigma = 1.3)), x), tolerance = 1e-12)
  expect_equal(mdens(marg("gauss", c(mu = 0.2, sigma = 1.3)), x), dnorm(x, 0.2, 1.3), tolerance = 1e-12)
  expect_equal(mdens(marg("st", c(mu = 0.2, sigma = 1.3, df = 4)), x), dt((x - 0.2) / 1.3, 4) / 1.3, tolerance = 1e-12)
  expect_equal(mdf(marg("st", c(mu = 0, sigma = 1, df = 1e6)), x), pnorm(x), tolerance = 1e-5)
  # Johnson SU tends to the Gaussian as delta grows with sigma = delta * s
  expect_equal(mdf(marg("jsu", c(mu = 0.2, sigma = 1e3 * 1.3, delta = 1e3, gamma = 0)), x),
               pnorm(x, 0.2, 1.3), tolerance = 1e-5)
  expect_equal(mdf(marg("jsu", c(mu = 0, sigma = 1, delta = 2, gamma = 0)), 0), 0.5)
  expect_error(marg("jsu", c(mu = 0, sigma = 1, delta = -1, gamma = 0)), "delta")
  # NIG tends to the Gaussian as shape grows
  expect_equal(mdf(marg("nig", c(mu = 0.2, sigma = 1.3, shape = 1e4, skew = 0)), x),
               pnorm(x, 0.2, 1.3), tolerance = 1e-4)
  expect_error(marg("nig", c(mu = 0, sigma = 1, shape = 0, skew = 0)), "shape")
  # gh contains the NIG (lambda = -1/2) and the hyperbolic distribution (lambda = 1)
  th <- c(mu = 0.2, sigma = 1.3, shape = 1.4, skew = 0.5)
  expect_equal(mdens(marg("gh", c(th[1:2], lambda = -0.5, th[3:4])), x), mdens(marg("nig", th), x), tolerance = 1e-9)
  expect_equal(mdens(marg("gh", c(th[1:2], lambda = 1, th[3:4])), x), mdens(marg("hyp", th), x), tolerance = 1e-9)
  expect_equal(mdf(marg("gh", c(th[1:2], lambda = -0.5, th[3:4])), x), mdf(marg("nig", th), x), tolerance = 1e-7)
  # the skew parameter moves mass to the right for gamma > 1
  expect_gt(mdf(marg("sst", c(mu = 0, sigma = 1, df = 6, gamma = 1.5)), 0, lower.tail = FALSE),
            mdf(marg("sst", c(mu = 0, sigma = 1, df = 6, gamma = 0.7)), 0, lower.tail = FALSE))
})

test_that("random draws follow the distribution and tssim/resid work", {
  set.seed(1)
  for (fam in names(mpars)) {
    m <- marg(fam, mpars[[fam]])
    x <- tssim(m, 4000)
    expect_length(x, 4000)
    r <- resid(m, x)
    expect_gt(stats::ks.test(r, "punif")$p.value, 0.001)
    expect_equal(resid(m, x, scale = "normal"), qnorm(r), tolerance = 1e-12)
    expect_equal(residuals(m, x), r)
  }
  xs <- ts(rnorm(30))
  expect_s3_class(resid(marg("gauss", c(mu = 0, sigma = 1)), xs), "ts")
})

test_that("the NIG random draws have the right mean and variance, and simulation is exact", {
  set.seed(2)
  th <- c(mu = 0.2, sigma = 1.3, shape = 1.5, skew = 0.4)
  x <- tssim(marg("nig", th), 2e5)
  expect_equal(mean(x), th[["mu"]] + th[["sigma"]] * th[["skew"]], tolerance = 0.01)
  expect_equal(var(x), th[["sigma"]]^2 * (1 + th[["skew"]]^2 / th[["shape"]]), tolerance = 0.03)
  expect_gt(stats::ks.test(mdf(marg("nig", th), x[1:5000]), "punif")$p.value, 0.001)
})

test_that("the NIG agrees with the ghyp package (gamma = sigma * skew)", {
  skip_if_not_installed("ghyp")
  th <- mpars$nig
  o <- ghyp::NIG(alpha.bar = th[["shape"]], mu = th[["mu"]], sigma = th[["sigma"]], gamma = th[["sigma"]] * th[["skew"]])
  m <- marg("nig", th)
  x <- seq(-9, 9, length.out = 31)
  p <- c(1e-8, 1e-3, 0.05, 0.5, 0.9, 0.999, 1 - 1e-6)
  expect_equal(mdens(m, x), ghyp::dghyp(x, o), tolerance = 1e-9)
  expect_equal(mdf(m, x), ghyp::pghyp(x, o), tolerance = 1e-7)
  expect_equal(mqf(m, p), ghyp::qghyp(p, o), tolerance = 1e-5)
})

test_that("hyp and gh agree with the ghyp package (gamma = sigma * skew)", {
  skip_if_not_installed("ghyp")
  p <- c(1e-8, 1e-3, 0.05, 0.5, 0.9, 0.999, 1 - 1e-6)
  x <- seq(-9, 9, length.out = 31)
  for (fam in c("hyp", "gh")) {
    th <- mpars[[fam]]
    lam <- if (fam == "gh") th[["lambda"]] else 1
    o <- ghyp::ghyp(lambda = lam, alpha.bar = th[["shape"]], mu = th[["mu"]], sigma = th[["sigma"]],
                    gamma = th[["sigma"]] * th[["skew"]])
    m <- marg(fam, th)
    expect_equal(mdens(m, x), ghyp::dghyp(x, o), tolerance = 1e-9, info = fam)
    expect_equal(mdf(m, x), ghyp::pghyp(x, o), tolerance = 1e-6, info = fam)
    expect_equal(mqf(m, p), ghyp::qghyp(p, o), tolerance = 1e-4, info = fam)
  }
})

test_that("margins agree with the legacy tscopula distribution functions", {
  skip_if_no_oracle()
  x <- seq(-3, 4, length.out = 25)
  p <- c(0.001, 0.05, 0.3, 0.5, 0.9, 0.995)
  th <- mpars$sst
  m <- marg("sst", th)
  expect_equal(mdens(m, x), tscopula::dsst(x, df = th[["df"]], gamma = th[["gamma"]], mu = th[["mu"]], sigma = th[["sigma"]]), tolerance = 1e-8)
  expect_equal(mdf(m, x), tscopula::psst(x, df = th[["df"]], gamma = th[["gamma"]], mu = th[["mu"]], sigma = th[["sigma"]]), tolerance = 1e-8)
  expect_equal(mqf(m, p), tscopula::qsst(p, df = th[["df"]], gamma = th[["gamma"]], mu = th[["mu"]], sigma = th[["sigma"]]), tolerance = 1e-8)
  th <- mpars$st
  m <- marg("st", th)
  expect_equal(mdens(m, x), tscopula::dst(x, df = th[["df"]], mu = th[["mu"]], sigma = th[["sigma"]]), tolerance = 1e-10)
  th <- mpars$laplace
  expect_equal(mdens(marg("laplace", th), x), tscopula::dlaplace(x, mu = th[["mu"]], scale = th[["sigma"]]), tolerance = 1e-10)
  expect_equal(mdf(marg("laplace", th), x), tscopula::plaplace(x, mu = th[["mu"]], scale = th[["sigma"]]), tolerance = 1e-10)
  th <- mpars$jsu
  mj <- marg("jsu", th)
  expect_equal(mdens(mj, x), tscopula::djsu(x, delta = th[["delta"]], mu = th[["mu"]], sigma = th[["sigma"]], gamma = th[["gamma"]]), tolerance = 1e-10)
  expect_equal(mdf(mj, x), tscopula::pjsu(x, delta = th[["delta"]], mu = th[["mu"]], sigma = th[["sigma"]], gamma = th[["gamma"]]), tolerance = 1e-10)
  expect_equal(mqf(mj, p), tscopula::qjsu(p, delta = th[["delta"]], mu = th[["mu"]], sigma = th[["sigma"]], gamma = th[["gamma"]]), tolerance = 1e-8)
  th <- mpars$slaplace
  expect_equal(mdens(marg("slaplace", th), x), tscopula::dslaplace(x, mu = th[["mu"]], scale = th[["sigma"]], gamma = th[["gamma"]]), tolerance = 1e-8)
  expect_equal(mdf(marg("slaplace", th), x), tscopula::pslaplace(x, mu = th[["mu"]], scale = th[["sigma"]], gamma = th[["gamma"]]), tolerance = 1e-8)
})
