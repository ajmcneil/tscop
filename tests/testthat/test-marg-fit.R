test_that("tsfit recovers the parameters of every family", {
  skip_on_cran()
  set.seed(2)
  truths <- list(
    gauss = c(mu = 0.3, sigma = 1.7),
    laplace = c(mu = -0.2, sigma = 0.8),
    slaplace = c(mu = 0.1, sigma = 1.2, gamma = 0.7),
    st = c(mu = 0.5, sigma = 0.9, df = 4.5),
    sst = c(mu = -0.3, sigma = 1.1, df = 6, gamma = 1.4),
    jsu = c(mu = 0.2, sigma = 1.5, delta = 2.2, gamma = -0.6)
  )
  for (fam in names(truths)) {
    x <- tssim(marg(fam, truths[[fam]]), 6000)
    f <- tsfit(x, marg(fam), se = "hessian")
    expect_s4_class(f, "margfit")
    expect_equal(f@convergence, 0L, info = fam)
    est <- coef(f)
    se <- f@se
    expect_true(all(abs(est - truths[[fam]]) < 4.5 * se), info = fam)
    expect_equal(f@se_method, "hessian")
    expect_true(all(is.finite(se) & se > 0))
    expect_s3_class(logLik(f), "logLik")
    expect_equal(attr(logLik(f), "df"), length(truths[[fam]]))
    expect_equal(as.numeric(logLik(f)), sum(mdens(f, x, log = TRUE)), tolerance = 1e-8)
  }
})

test_that("fixed parameters are respected and the likelihood is no larger than the free fit", {
  set.seed(3)
  x <- tssim(marg("st", c(mu = 0, sigma = 1, df = 5)), 3000)
  ff <- tsfit(x, marg("st"))
  fx <- tsfit(x, marg("st", c(mu = 0), fix = "mu"))
  expect_equal(fx@model@pars[["mu"]], 0)
  expect_equal(names(coef(fx)), c("sigma", "df"))
  expect_equal(attr(logLik(fx), "df"), 2L)
  expect_lte(as.numeric(logLik(fx)), as.numeric(logLik(ff)) + 1e-6)
  # fixing the skewing parameter of the sst at one gives the st fit
  fs <- tsfit(x, marg("sst", c(gamma = 1), fix = "gamma"))
  expect_equal(as.numeric(logLik(fs)), as.numeric(logLik(ff)), tolerance = 1e-5)
})

test_that("init, bootstrap standard errors and input checks", {
  skip_on_cran()
  set.seed(4)
  x <- tssim(marg("sst", c(mu = 0, sigma = 1, df = 6, gamma = 1.2)), 1500)
  f <- tsfit(x, marg("sst"), se = "hessian")
  g <- tsfit(x, marg("sst"), init = f@model)
  expect_equal(as.numeric(logLik(g)), as.numeric(logLik(f)), tolerance = 1e-6)
  fb <- tsfit(x, marg("sst"), se = "bootstrap", B = 30)
  expect_equal(fb@se_method, "bootstrap")
  expect_equal(unname(fb@se), unname(f@se), tolerance = 0.5)
  expect_error(tsfit(x, marg("sst"), se = "foo"), "se")
  expect_error(tsfit(x[1:5], marg("sst")), "at least 10")
  expect_error(tsfit("a", marg("sst")), "numeric")
  expect_equal(unname(coef(tsfit(c(x, NA, Inf), marg("gauss")))), unname(coef(tsfit(x, marg("gauss")))))
  expect_output(show(f), "standard errors: hessian")
  expect_equal(length(tssim(f, 50)), 50)
  expect_equal(length(resid(f)), length(x))
  expect_gt(stats::ks.test(resid(f), "punif")$p.value, 0.001)
})
