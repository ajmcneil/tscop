test_that("plot draws the pair copulas of an sdvinecop model", {
  m <- sdvinecop(arfactcop(lambda = 0.9, w = 0.5), tail = pc("ivt", base = "joe"), vt = c(0.45, 0.65), vttype = "v2p",
                 vtshape = c(0.7, 1.4), maxlag = 8, tautol = 1e-3)
  pdf(NULL); on.exit(dev.off())
  r <- plot(m)
  expect_equal(r$lags, 1:4)
  expect_length(r$tau, 4)
  expect_equal(r$tau, ktau(m)[1:4], tolerance = 1e-8)
  expect_true(all(diff(r$levels) > 0))
  expect_silent(plot(m, lags = c(1, 3), strips = TRUE))
  expect_error(plot(m, lags = 99), "effective lag")
  expect_error(plot(m, points = 1), "'data' is needed")
})

test_that("the data pairs and the two-model layout work", {
  set.seed(5)
  m <- sdvinecop(arfactcop(lambda = 0.9, w = 0.5), tail = pc("ivt", base = "joe"), vt = c(0.45, 0.65), maxlag = 8, tautol = 1e-3)
  mn <- sdvinecop(arfactcop(lambda = 0.9, w = 0.5), tail = pc("ivt", base = "joe"), vt = c(0.45, 0.65), vttype = "v2p",
                  vtshape = c(0.7, 1.4), maxlag = 8, tautol = 1e-3)
  u <- tssim(m, 400)
  pdf(NULL); on.exit(dev.off())
  expect_silent(plot(m, lags = 1:2, data = u))
  expect_silent(plot(m, mn, lags = 1, data = u, strips = TRUE, labels = c("linear", "v2p")))
  # the level pairs at lag 1 are the consecutive observations; the plot uses the levels pooled over both models
  pp <- dv_level_tau_r(dv_resolve(m)[1:2], u, keep = 2)$pairs
  expect_equal(unname(pp[[1]]), cbind(u[-400], u[-1]))
  r <- plot(m, mn, lags = 1:2)
  expect_equal(dim(r$tau), c(2L, 2L))
  expect_error(plot(m, data = c(0.5, 1.2, 0.3, 0.4, 0.2, 0.9, 0.8), points = 1), "strictly inside")
})

test_that("plot methods for fits use the data of the fit", {
  skip_on_cran()
  set.seed(9)
  tru <- tsmod(marg("st", c(mu = 0, sigma = 1, df = 6)),
               sdvinecop(arfactcop(lambda = 0.9, w = 0.5), tail = pc("ivt", base = "joe"), vt = c(0.45, 0.6), maxlag = 8, tautol = 1e-3))
  x <- tssim(tru, 500)
  ft <- tsfit(x, tsmod("st", sdvinecop(arfactcop(order = 1), tail = pc("ivt", base = "joe"), maxlag = 8, tautol = 1e-3)))
  pdf(NULL); on.exit(dev.off())
  expect_silent(plot(ft, lags = 1:2, data = TRUE, strips = TRUE))
  fj <- tsfit(x, tsmod(marg("st"), sdvinecop(arfactcop(order = 1), tail = pc("ivt", base = "joe"), maxlag = 8, tautol = 1e-3)), method = "joint")
  expect_silent(plot(fj, lags = 1, data = TRUE))
})
