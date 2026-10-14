test_that("joint ML with Gaussian margin and ARMA copula is the exact Gaussian ARMA likelihood", {
  set.seed(1)
  x <- 1 + 2 * as.numeric(arima.sim(list(ar = 0.6, ma = 0.3), 600))
  a <- arima(x, order = c(1, 0, 1), method = "ML")
  f <- tsfit(x, tsmod("gauss", armacop(order = c(1, 1))), method = "joint")
  expect_equal(as.numeric(logLik(f)), as.numeric(a$loglik), tolerance = 1e-5)
  expect_equal(unname(f@estimate["copula.ar1"]), unname(a$coef["ar1"]), tolerance = 1e-3)
  expect_equal(unname(f@estimate["copula.ma1"]), unname(a$coef["ma1"]), tolerance = 1e-3)
  expect_equal(unname(f@estimate["margin.mu"]), unname(a$coef["intercept"]), tolerance = 1e-3)
  expect_equal(f@method, "joint")
})

test_that("joint ML never has a lower likelihood than the two-stage fit and recovers the model", {
  set.seed(2)
  tru <- tsmod(marg("sst", c(mu = 0, sigma = 1, df = 6, gamma = 1.2)), arfactcop(lambda = 0.9, w = 0.6))
  x <- tssim(tru, 1500)
  f2 <- tsfit(x, tsmod("sst", arfactcop(order = 1)))
  fj <- tsfit(x, tsmod("sst", arfactcop(order = 1)), method = "joint")
  expect_gte(as.numeric(logLik(fj)), as.numeric(logLik(f2)) - 1e-6)
  expect_equal(names(coef(fj)), names(coef(f2)))
  expect_equal(unname(coef(fj)), unname(coef(f2)), tolerance = 0.15)
  # the reported likelihood is margin + copula likelihood at the joint estimate
  u <- tm_pits(fj@model@margin, x)
  cl <- -cop_obj(fj@model@copula, u)$negll(cop_obj(fj@model@copula, u)$z0, u)
  expect_equal(as.numeric(logLik(fj)), sum(mdens(fj@model@margin, x, log = TRUE)) + cl, tolerance = 1e-6)
  expect_equal(fj@npar, length(coef(fj)))
})

test_that("cop_obj reproduces the log-likelihood of every copula class", {
  set.seed(3)
  x <- tssim(tsmod(marg("st", c(mu = 0, sigma = 1, df = 5)), armacop(ar = 0.6, ma = 0.3)), 500)
  u <- rank(x) / (length(x) + 1)
  fits <- list(
    tsfit(u, armacop(order = c(1, 1))),
    tsfit(u, arfactcop(order = 1)),
    tsfit(u, sdvinecop(armacop(order = c(1, 1)), pairs = list(pc("gumbel")), tail = "gumbel")),
    tsfit(u, vtscop(arfactcop(order = 1), udp::v2p(0.45, 0.8)))
  )
  for (f in fits) {
    o <- cop_obj(f@model, u)
    expect_equal(-o$negll(o$z0, u), as.numeric(logLik(f)), tolerance = 1e-6)
    expect_equal(unname(o$coef(o$z0)), unname(coef(f)), tolerance = 1e-6)
    expect_equal(coef(o$build(o$z0)), coef(f@model))
  }
})

test_that("joint ML works for d-vine and v-transform copulas and gives Hessian standard errors", {
  skip_on_cran()
  set.seed(4)
  tru <- tsmod(marg("st", c(mu = 0, sigma = 1, df = 5)),
               sdvinecop(armacop(ar = 0.7, ma = 0.2), pairs = list(pc("gumbel")), tail = "gumbel"))
  x <- tssim(tru, 800)
  spec <- tsmod("st", sdvinecop(armacop(order = c(1, 1)), pairs = list(pc("gumbel")), tail = "gumbel"))
  f2 <- tsfit(x, spec)
  fj <- tsfit(x, spec, method = "joint", se = "hessian")
  expect_gte(as.numeric(logLik(fj)), as.numeric(logLik(f2)) - 1e-6)
  expect_equal(fj@se_method, "hessian")
  expect_true(all(fj@se > 0))
  expect_equal(fj@model@copula@K, f2@model@copula@K)
  expect_output(show(fj), "method: joint")
})

test_that("joint ML is refused for an empirical margin and bootstrap SEs work", {
  skip_on_cran()
  set.seed(5)
  x <- tssim(tsmod(marg("st", c(mu = 0, sigma = 1, df = 5)), arfactcop(lambda = 0.9, w = 0.6)), 400)
  expect_error(tsfit(x, tsmod(emarg(x), arfactcop(order = 1)), method = "joint"), "parametric margin")
  expect_error(tsfit(x, tsmod("st", arfactcop(order = 1)), se = "hessian"), "joint")
  fb <- tsfit(x, tsmod("st", arfactcop(order = 1)), method = "joint", se = "bootstrap", B = 10)
  expect_equal(fb@se_method, "bootstrap")
  expect_equal(dim(fb@boot), c(10L, 5L))
})

test_that("the exact joint gradient equals finite differences (d-vine, IVT and v-transform copulas)", {
  set.seed(12)
  kp <- arfactcop(lambda = 0.9, w = 0.5)
  cases <- list(
    list(tsmod(marg("st", c(mu = 0, sigma = 1, df = 6)),
               sdvinecop(kp, pairs = list(pc("gumbel")), tail = pc("joe"), maxlag = 12, tautol = 1e-4)),
         tsmod("st", sdvinecop(arfactcop(order = 1), pairs = list(pc("gumbel")), tail = pc("joe"), maxlag = 12, tautol = 1e-4))),
    list(tsmod(marg("sst", c(mu = 0, sigma = 1, df = 6, gamma = 1.2)),
               sdvinecop(kp, tail = pc("ivt", base = "joe"), vt = c(0.45, 0.55), maxlag = 12, tautol = 1e-4)),
         tsmod("sst", sdvinecop(arfactcop(order = 1), tail = pc("ivt", base = "joe"), maxlag = 12, tautol = 1e-4))),
    list(tsmod(marg("st", c(mu = 0, sigma = 1, df = 8)),
               vtscop(sdvinecop(kp, pairs = list(pc("joe")), tail = pc("joe"), maxlag = 12, tautol = 1e-4), udp::v2p(0.45, 0.8))),
         tsmod("st", vtscop(sdvinecop(arfactcop(order = 1), pairs = list(pc("joe")), tail = pc("joe"), maxlag = 12, tautol = 1e-4),
                            udp::v2p(0.5, 1))))
  )
  for (cs in cases) {
    x <- tssim(cs[[1]], 600)
    f2 <- tsfit(x, cs[[2]])
    jt <- tm_joint_setup(f2@model@margin, f2@model@copula, x)
    th <- jt$theta0 + rnorm(length(jt$theta0), 0, 0.02)
    g <- jt$gradient(th, x)
    expect_false(is.null(g))
    fd <- vapply(seq_along(th), function(k) {
      h <- 1e-6
      tp <- th; tp[k] <- tp[k] + h
      tm <- th; tm[k] <- tm[k] - h
      (jt$negll(tp, x) - jt$negll(tm, x)) / (2 * h)
    }, 1)
    expect_lt(max(abs(g - fd) / pmax(abs(fd), 1)), 1e-4)
  }
})

test_that("joint fits with the exact gradient agree with fits by differences", {
  skip_on_cran()
  set.seed(13)
  kp <- arfactcop(lambda = 0.9, w = 0.5)
  tru <- tsmod(marg("st", c(mu = 0, sigma = 1, df = 6)),
               sdvinecop(kp, pairs = list(pc("gumbel")), tail = pc("joe"), maxlag = 10, tautol = 1e-4))
  spec <- tsmod("st", sdvinecop(arfactcop(order = 1), pairs = list(pc("gumbel")), tail = pc("joe"), maxlag = 10, tautol = 1e-4))
  x <- tssim(tru, 700)
  f1 <- tsfit(x, spec, method = "joint")
  old <- options(tscop.cppgrad = FALSE)
  on.exit(options(old))
  f2 <- tsfit(x, spec, method = "joint")
  expect_equal(unname(coef(f1)), unname(coef(f2)), tolerance = 1e-3)
  expect_equal(as.numeric(logLik(f1)), as.numeric(logLik(f2)), tolerance = 1e-6)
})
