make_fc <- function(p = c(0.01, 0.05)) {
  set.seed(7)
  m <- tsmod(marg("st", c(mu = 0, sigma = 1, df = 6)), armacop(ar = 0.5))
  x <- tssim(m, 400)
  g <- tsmod(marg("gauss", c(mu = 0, sigma = 1.3)), armacop(ar = 0.5))
  list(good = tsroll(m, x, start = 301, p = p), bad = tsroll(g, x, start = 301, p = p))
}

test_that("tsroll with p = NULL has no risk measures and the same forecasts", {
  fc <- make_fc(); f0 <- make_fc(NULL)
  expect_equal(f0$good@p, numeric(0))
  expect_equal(dim(f0$good@var), c(100L, 0L))
  expect_equal(dim(f0$good@es), c(100L, 0L))
  expect_equal(f0$good@quantile, fc$good@quantile)
  expect_equal(f0$good@logdens, fc$good@logdens)
  expect_output(show(f0$good), "100 forecasts")
  sc <- with(f0, tsscore(bad, good))
  expect_named(sc$summary, c("forecast", "aqs", "neglogscore"))
  expect_equal(length(sc$var), 0L)
  expect_error(tsscore(f0$bad, p = 0.05), "contain no value-at-risk")
})

test_that("as.data.frame.tsscore is a tidy version of the summary and the tests", {
  fc <- make_fc()
  sc <- with(fc, tsscore(bad, good))
  d <- as.data.frame(sc)
  expect_named(d, c("forecast", "loss", "score", "change", "change_pct", "statistic", "p.value"))
  expect_equal(attr(d, "reference"), "bad")
  expect_setequal(unique(d$loss), c("aqs", "neglogscore", "fz0_0.01", "fz0_0.05"))
  a <- d[d$loss == "aqs", ]
  expect_equal(a$score, sc$summary$aqs)
  expect_equal(a$change, c(0, sc$summary$aqs[2] - sc$summary$aqs[1]))
  expect_equal(a$change_pct[2], 100 * (sc$summary$aqs[2] / sc$summary$aqs[1] - 1))
  dm <- sc$dm[sc$dm$loss == "aqs", ]
  expect_equal(a$statistic, c(NA, dm$statistic))
  expect_equal(a$p.value, c(NA, dm$p.value))
  # a smaller loss and a positive statistic go together
  expect_lt(a$change[2], 0); expect_gt(a$statistic[2], 0)
  # the log score can be negative: no percentage when the reference score is not positive
  l <- d[d$loss == "neglogscore", ]
  expect_equal(l$score, sc$summary$neglogscore)
})

test_that("tstable arranges one loss for one or several set-ups", {
  fc <- make_fc()
  sc <- with(fc, tsscore(bad, good))
  t1 <- tstable(sc, scale = 1000)
  expect_s3_class(t1, "tstable")
  expect_equal(t1$forecast, c("bad", "good"))
  expect_equal(t1$score, 1000 * sc$summary$aqs)
  expect_equal(t1$change[2], 100 * (sc$summary$aqs[2] / sc$summary$aqs[1] - 1))
  expect_true(is.na(t1$p.value[1]))
  expect_output(print(t1), "Loss: aqs \\(scores multiplied by 1000\\)")
  expect_output(print(t1), "percentage change relative to 'bad'")
  # the difference for the log score
  t2 <- tstable(sc, loss = "neglogscore")
  expect_equal(attr(t2, "relative"), "difference")
  expect_equal(t2$change[2], sc$summary$neglogscore[2] - sc$summary$neglogscore[1])
  expect_equal(tstable(sc, relative = "difference")$change[2], sc$summary$aqs[2] - sc$summary$aqs[1])
  # several set-ups side by side
  t3 <- tstable(list(fixed = sc, other = sc))
  expect_named(t3, c("forecast", "score.fixed", "change.fixed", "p.value.fixed", "score.other", "change.other", "p.value.other"))
  expect_output(print(t3), "score.fixed")
  # errors
  expect_error(tstable(sc, loss = "nothing"), "no loss 'nothing'")
  expect_error(tstable(list(sc, sc)), "must be named")
  expect_error(tstable(list(a = sc, b = with(fc, tsscore(good, bad)))), "must be the same")
  expect_error(tstable(1), "tsscore object")
})
