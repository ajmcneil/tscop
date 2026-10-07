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
