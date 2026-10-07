# Polynomial helpers for (seasonal) ARMA copula processes. Conventions follow
# stats::arima: x_t = sum ar_i x_{t-i} + e_t + sum ma_j e_{t-j}, i.e. the AR
# polynomial is 1 - sum ar_i B^i and the MA polynomial is 1 + sum ma_j B^j.

# Product of two polynomials given by coefficient vectors (constant term first).
poly_mult <- function(a, b) {
  out <- numeric(length(a) + length(b) - 1L)
  for (i in seq_along(a)) {
    idx <- i - 1L + seq_along(b)
    out[idx] <- out[idx] + a[i] * b
  }
  out
}

# Coefficients of a seasonal polynomial in B, constant term first: the
# seasonal coefficients 'coef' sit at lags period, 2 * period, ...
seasonal_poly <- function(coef, period, sign) {
  out <- numeric(length(coef) * period + 1L)
  out[1L] <- 1
  out[period * seq_along(coef) + 1L] <- sign * coef
  out
}

# Expand multiplicative seasonal AR/MA coefficients into the equivalent plain
# ARMA coefficients. Returns list(ar, ma), each possibly numeric(0).
arma_expand <- function(ar = numeric(0), ma = numeric(0),
                        sar = numeric(0), sma = numeric(0), period = NA) {
  if ((length(sar) > 0L || length(sma) > 0L) &&
      (is.na(period) || period < 2)) {
    stop("'period' must be at least 2 when seasonal terms are present")
  }
  arpoly <- c(1, -ar)
  mapoly <- c(1, ma)
  if (length(sar) > 0L) {
    arpoly <- poly_mult(arpoly, seasonal_poly(sar, period, -1))
  }
  if (length(sma) > 0L) {
    mapoly <- poly_mult(mapoly, seasonal_poly(sma, period, 1))
  }
  list(ar = -arpoly[-1L], ma = mapoly[-1L])
}

# Partial autocorrelations of an AR polynomial (backward Durbin-Levinson).
# Stationary if and only if all |pacf| < 1; an exact unit PACF value stops the
# recursion and is returned as is.
ar_to_pacf <- function(ar) {
  p <- length(ar)
  pacf <- numeric(p)
  phi <- ar
  for (k in rev(seq_len(p))) {
    pacf[k] <- phi[k]
    if (k == 1L) break
    if (abs(phi[k]) >= 1) {
      pacf[seq_len(k - 1L)] <- NA_real_
      break
    }
    phi <- (phi[seq_len(k - 1L)] + phi[k] * rev(phi[seq_len(k - 1L)])) /
      (1 - phi[k]^2)
  }
  pacf
}

# AR coefficients from partial autocorrelations (forward Durbin-Levinson).
pacf_to_ar <- function(pacf) {
  p <- length(pacf)
  phi <- numeric(0)
  for (k in seq_len(p)) {
    phi <- c(phi - pacf[k] * rev(phi), pacf[k])
  }
  phi
}

# TRUE if the AR polynomial 1 - sum ar_i B^i has all roots outside the unit
# circle (a causal, stationary process). The empty polynomial is stationary.
arma_is_stationary <- function(ar) {
  if (length(ar) == 0L) return(TRUE)
  pacf <- ar_to_pacf(ar)
  !anyNA(pacf) && all(abs(pacf) < 1)
}

# TRUE if the MA polynomial 1 + sum ma_j B^j is invertible; this is the
# stationarity condition for the polynomial with coefficients -ma.
arma_is_invertible <- function(ma) {
  arma_is_stationary(-ma)
}

# Unconstrained <-> admissible maps. An unconstrained vector is mapped by
# tanh to partial autocorrelations in (-1, 1) and then to AR coefficients, so
# any real vector gives a stationary polynomial. The same map is used for MA
# coefficients, via -ma.
free_to_ar <- function(z) {
  pacf_to_ar(tanh(z))
}

ar_to_free <- function(ar) {
  atanh(ar_to_pacf(ar))
}

free_to_ma <- function(z) {
  -free_to_ar(z)
}

ma_to_free <- function(ma) {
  ar_to_free(-ma)
}
