test_that("emarg validates its arguments", {
  set.seed(1)
  x <- rt(200, 5)
  e <- emarg(x)
  expect_s4_class(e, "emarg")
  expect_equal(e@data, x)
  expect_equal(e@sorted, sort(x))
  expect_equal(emarg(c(x, NA, Inf))@data, x)
  expect_error(emarg(x[1:3]), "at least 5")
  expect_error(emarg("a"), "numeric")
  expect_error(emarg(x, qtype = 10), "qtype")
  expect_error(emarg(x, bw = -1), "bw")
  expect_error(emarg(x, type = "foo"))
  expect_equal(emarg(x, bw = 0.3)@bw, 0.3)
  expect_equal(emarg(x, bw = "nrd0")@bw, bw.nrd0(x))
  expect_output(show(e), "emarg")
})

test_that("the scaled edf gives rank / (n + 1) at the data and stays inside (0, 1)", {
  set.seed(2)
  x <- rnorm(150)
  e <- emarg(x)
  n <- length(x)
  expect_equal(mdf(e, x), rank(x) / (n + 1))
  expect_equal(resid(e), rank(x) / (n + 1))
  expect_equal(mdf(e, c(-100, 100)), c(0.5 / (n + 1), (n + 0.5) / (n + 1)))
  expect_equal(mdf(e, x, lower.tail = FALSE), 1 - rank(x) / (n + 1))
  ep <- emarg(x, proper = TRUE)
  expect_equal(mdf(ep, c(sort(x), 100)), c(seq_len(n) / n, 1))
  expect_equal(resid(e, scale = "normal"), qnorm(rank(x) / (n + 1)))
  expect_equal(residuals(e), resid(e))
})

test_that("quantile functions of the edf type are the empirical quantiles", {
  set.seed(3)
  x <- rnorm(100)
  p <- c(0, 0.01, 0.3, 0.5, 0.97, 1)
  for (t in c(1, 4, 7, 8)) {
    expect_equal(mqf(emarg(x, qtype = t), p), as.numeric(quantile(x, p, type = t)))
  }
  expect_equal(mqf(emarg(x), c(0.1, 0.9), lower.tail = FALSE), as.numeric(quantile(x, c(0.9, 0.1), type = 7)))
  expect_error(mqf(emarg(x), 1.5), "\\[0, 1\\]")
})

test_that("the kernel margin is a smooth distribution: density, df and quantile function agree", {
  set.seed(4)
  x <- rt(200, 4)
  e <- emarg(x, "kernel")
  g <- seq(min(x) - 5 * e@bw, max(x) + 5 * e@bw, length.out = 41)
  expect_equal(integrate(function(v) mdens(e, v), -Inf, Inf, rel.tol = 1e-8)$value, 1, tolerance = 1e-6)
  h <- 1e-5
  expect_equal((mdf(e, g + h) - mdf(e, g - h)) / (2 * h), mdens(e, g), tolerance = 1e-5)
  expect_equal(mdf(e, g) + mdf(e, g, lower.tail = FALSE), rep(1, length(g)), tolerance = 1e-12)
  p <- c(1e-6, 0.01, 0.3, 0.5, 0.9, 0.999)
  expect_equal(mdf(e, mqf(e, p)), p, tolerance = 1e-9)
  expect_equal(mqf(e, p, lower.tail = FALSE), mqf(e, 1 - p), tolerance = 1e-8)
  expect_equal(mqf(e, c(0, 1)), c(-Inf, Inf))
  expect_equal(mdens(e, 0, log = TRUE), log(mdens(e, 0)))
  # the density of the edf type is the same kernel estimate
  expect_equal(mdens(emarg(x, "edf"), g), mdens(e, g))
  # it matches stats::density
  d <- density(x, bw = e@bw, kernel = "gaussian", n = 1024)
  expect_equal(mdens(e, d$x[c(100, 300, 700)]), d$y[c(100, 300, 700)], tolerance = 5e-3)
})

test_that("resampling and the time series template", {
  set.seed(5)
  x <- rnorm(80)
  e <- emarg(x)
  s <- tssim(e, 200)
  expect_length(s, 200)
  expect_true(all(s %in% x))
  expect_length(tssim(e), 80)
  sk <- tssim(emarg(x, "kernel"), 300)
  expect_false(any(sk %in% x))
  xs <- ts(x, start = 2000)
  es <- emarg(xs)
  expect_s3_class(resid(es), "ts")
  expect_equal(as.numeric(resid(es)), rank(x) / (81))
})

test_that("a copula process can be fitted to the pseudo-observations of an empirical margin", {
  skip_on_cran()
  set.seed(6)
  u <- tssim(arfactcop(lambda = 0.8, w = 0.6), 800)
  x <- qt(u, 5)                                     # observations with a t margin and this dependence
  u_hat <- resid(emarg(x))                          # pseudo-observations
  expect_gt(cor(u, u_hat, method = "spearman"), 0.999)
  f <- tsfit(u_hat, arfactcop(order = 1), pseudo = FALSE)
  expect_equal(unname(coef(f)), c(0.8, 0.6), tolerance = 0.15)
})
