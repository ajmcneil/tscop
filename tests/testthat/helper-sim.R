# Test-only simulator of the unit-variance Gaussian ARMA normal scores, via
# arima.sim with a long burn-in (tssim will provide an exact simulator).
sim_arma_scores <- function(n, ar = numeric(0), ma = numeric(0), seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  g0 <- arma_acvf(ar, ma, 0)
  model <- list()
  if (length(ar) > 0) model$ar <- ar
  if (length(ma) > 0) model$ma <- ma
  as.numeric(arima.sim(model, n = n, sd = 1 / sqrt(g0), n.start = 2000))
}

# Test-only simulator of the unit-variance Gaussian factor process (sum of
# independent stationary AR(1)s plus noise): normal scores.
sim_arfact_scores <- function(n, lambda, w, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  w0 <- max(1 - sum(w), 0)
  x <- sqrt(w0) * rnorm(n)
  for (i in seq_along(lambda)) {
    e <- sqrt(1 - lambda[i]^2) * rnorm(n)
    e[1] <- rnorm(1)
    y <- as.numeric(stats::filter(e, lambda[i], method = "recursive"))
    x <- x + sqrt(w[i]) * y
  }
  x
}
