# stats::ARMAacf rejects the empty model; white noise has acf (1, 0, 0, ...).
acf_ref <- function(ar, ma, lag.max) {
  if (length(ar) == 0 && length(ma) == 0) return(c(1, numeric(lag.max)))
  as.numeric(ARMAacf(ar = ar, ma = ma, lag.max = lag.max))
}

# Conditional mean and sd of x_{n+h} given the observed x by Gaussian
# conditioning with the Toeplitz correlation matrix (NA positions skipped).
cond_moments <- function(x, acf, h) {
  n <- length(x)
  obs <- which(!is.na(x))
  Rfull <- toeplitz(acf)
  S <- Rfull[obs, obs, drop = FALSE]
  r <- acf[abs((n + h) - obs) + 1]
  w <- solve(S, r)
  c(mean = sum(w * x[obs]), sd = sqrt(1 - sum(w * r)))
}
