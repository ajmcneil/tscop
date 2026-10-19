test_that("quantile scores, AQS and VaR tests are computed correctly", {
  set.seed(3)
  n <- 2000
  y <- rnorm(n)
  lev <- (1:19) / 20
  Q <- matrix(qnorm(lev), n, length(lev), byrow = TRUE)
  Q2 <- Q + 0.5
  s <- tsscore(Q, shifted = Q2, observed = y, levels = lev, p = 0.05)
  pin <- function(q, tau) mean((as.numeric(y < q) - tau) * (q - y))
  prof <- vapply(seq_along(lev), function(j) pin(Q[1L, j], lev[j]), 0)
  expect_equal(unname(s$profile[, 1L]), prof, tolerance = 1e-12)
  expect_equal(s$summary$aqs[1L], mean(prof), tolerance = 1e-12)
  # the correct forecast beats the shifted one and the DM test detects it
  expect_lt(s$summary$aqs[1L], s$summary$aqs[2L])
  expect_lt(subset(s$dm, loss == "aqs")$statistic, -3)
  # Kupiec for the correct forecast does not reject; for the shifted one it does
  expect_gt(s$var[[1L]]$kupiec.p[1L], 0.01)
  expect_equal(s$var[[1L]]$breaches[1L], sum(y < qnorm(0.05)))
  expect_lt(s$var[[1L]]$kupiec.p[2L], 1e-3)
  # weights
  s2 <- tsscore(Q, observed = y, levels = lev, weights = function(l) l)
  expect_equal(s2$summary$aqs, sum(prof * lev / sum(lev)), tolerance = 1e-12)
  expect_output(print(s), "Diebold-Mariano")
})

test_that("Kupiec statistic matches the closed form", {
  I <- c(rep(0L, 95), rep(1L, 5)); p <- 0.01; n <- 100; x <- 5
  lr <- -2 * (x * log(p) + (n - x) * log(1 - p) - x * log(x / n) - (n - x) * log(1 - x / n))
  expect_equal(unname(ts_var_tests(I, p)["kupiec"]), lr, tolerance = 1e-12)
  # clustered breaches are detected by the independence test
  J <- c(rep(0L, 90), rep(1L, 10))
  expect_lt(ts_var_tests(J, 0.1)[["independence.p"]], 0.05)
})

test_that("ES statistics are near zero for correct forecasts and flag underestimation", {
  skip_on_cran()
  set.seed(4)
  m <- tsmod(marg("gauss", c(mu = 0, sigma = 1)),
             sdvinecop(arfactcop(lambda = 0.9, w = 0.3), tail = pc("ivt", base = "joe"), vt = c(0.5, 0.5), maxlag = 6, tautol = 1e-3))
  x <- tssim(m, 1200)
  fc <- tsroll(m, x, start = 201, p = 0.05)
  s <- tsscore(fc)
  e <- s$es[[1L]]
  expect_lt(abs(e$z2), 0.3)
  expect_gt(e$z2.p, 0.01)
  # shrinking the ES towards zero (underestimating risk) gives a positive Z2 and a worse FZ0
  fc2 <- fc; fc2@es <- fc@es * 0.7
  s2 <- tsscore(fc, bad = fc2)
  expect_gt(s2$es[[1L]]$z2[2L], s2$es[[1L]]$z2[1L])
  expect_gt(s2$es[[1L]]$fz0[2L], s2$es[[1L]]$fz0[1L])
})

test_that("tsscore validates its input", {
  y <- rnorm(10); Q <- matrix(0, 10, 3)
  expect_error(tsscore(Q), "observed")
  expect_error(tsscore(Q, observed = y), "levels")
  expect_error(tsscore(Q, observed = y, levels = c(0.1, 0.5)), "columns")
  expect_error(tsscore(Q, other = Q[1:5, ], observed = y, levels = c(.1, .5, .9)), "one row")
  expect_error(tsscore(Q, observed = y, levels = c(.1, .5, .9), weights = c(1, 2)), "weights")
})

test_that("the plot method runs", {
  set.seed(5); y <- rnorm(50); lev <- (1:9) / 10
  Q <- matrix(qnorm(lev), 50, 9, byrow = TRUE)
  s <- tsscore(Q, b = Q + 0.1, observed = y, levels = lev)
  pdf(NULL); on.exit(dev.off())
  expect_no_error(plot(s)); expect_no_error(plot(s, relative = TRUE))
})

test_that("forecasts are labelled by their names or argument names", {
  y <- c(0.1, -0.2, 0.3, 0); lev <- c(0.25, 0.75)
  Q1 <- matrix(c(-0.5, 0.5), 4, 2, byrow = TRUE); Q2 <- Q1 * 2
  s <- tsscore(Q1, wide = Q2, Q2, observed = y, levels = lev)
  expect_equal(s$summary$forecast, c("Q1", "wide", "Q2"))
  s2 <- tsscore(Q1 * 1, Q2 * 1, observed = y, levels = lev)
  expect_equal(s2$summary$forecast, c("forecast1", "forecast2"))
})
