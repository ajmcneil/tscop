# Autocovariances, state-space filter and likelihood for the zero-mean,
# unit-variance Gaussian ARMA process underlying an ARMA copula. All functions
# take plain (already expanded) ARMA coefficients in the stats::arima sign
# convention.

# Autocovariances gamma_0, ..., gamma_lag.max of the ARMA process driven by
# unit-variance innovations. Solves the (p + 1) x (p + 1) linear system for
# gamma_0, ..., gamma_p implied by the ARMA equations, then extends by
# recursion. No Kronecker products and no dependency on other packages.
arma_acvf <- function(ar = numeric(0), ma = numeric(0), lag.max = 0L) {
  p <- length(ar)
  q <- length(ma)
  # psi weights of the MA(infinity) form, psi_0, ..., psi_q
  psi <- numeric(q + 1L)
  psi[1L] <- 1
  for (j in seq_len(q)) {
    s <- ma[j]
    for (i in seq_len(min(j, p))) s <- s + ar[i] * psi[j - i + 1L]
    psi[j + 1L] <- s
  }
  # c_k = sum_{j = k}^q ma_j psi_{j - k}, with ma_0 = 1
  th <- c(1, ma)
  ck <- function(k) {
    if (k > q) return(0)
    sum(th[(k + 1L):(q + 1L)] * psi[seq_len(q + 1L - k)])
  }
  # system for gamma_0, ..., gamma_p
  A <- diag(p + 1L)
  for (k in 0:p) {
    for (i in seq_len(p)) {
      col <- abs(k - i) + 1L
      A[k + 1L, col] <- A[k + 1L, col] - ar[i]
    }
  }
  rhs <- vapply(0:p, ck, numeric(1))
  gam <- solve(A, rhs)
  nout <- max(lag.max, p) + 1L
  out <- numeric(nout)
  out[seq_len(p + 1L)] <- gam
  for (k in seq_len(nout - 1L)) {
    if (k <= p) next
    out[k + 1L] <- sum(ar * out[k - seq_len(p) + 1L]) + ck(k)
  }
  out[seq_len(lag.max + 1L)]
}

# Stationary covariance P of the state alpha_t = Tm alpha_{t-1} + R e_t, the
# solution of P = Tm P Tm' + RR, by the doubling algorithm: O(m^3 log) work
# rather than an m^2 x m^2 linear solve.
arma_stationary_cov <- function(Tm, RR, maxit = 200L) {
  P <- RR
  A <- Tm
  for (it in seq_len(maxit)) {
    dP <- A %*% P %*% t(A)
    P <- P + dP
    A <- A %*% A
    if (max(abs(dP)) <= 1e-16 * max(abs(P))) break
  }
  P
}

# State-space form (Harvey) of the zero-mean, unit-variance Gaussian ARMA
# process: alpha_t = Tm alpha_{t-1} + R e_t with e_t ~ N(0, 1), x_t = alpha_t[1].
# R carries the factor 1 / sqrt(gamma_0) that fixes Var(x_t) = 1. P0 is the
# stationary state covariance.
arma_state_space <- function(ar, ma) {
  p <- length(ar)
  q <- length(ma)
  m <- max(p, q + 1L)
  phi <- c(ar, numeric(m - p))
  theta <- c(1, ma, numeric(m - 1L - q))
  g0 <- arma_acvf(ar, ma, 0L)
  R <- theta / sqrt(g0)
  Tm <- matrix(0, m, m)
  Tm[, 1L] <- phi
  if (m > 1L) Tm[cbind(seq_len(m - 1L), seq_len(m - 1L) + 1L)] <- 1
  RR <- tcrossprod(R)
  list(m = m, phi = phi, Tm = Tm, R = R, RR = RR,
       P0 = arma_stationary_cov(Tm, RR))
}

# Exact simulation of n values of the process: the first state is drawn from
# its stationary distribution, so no burn-in is needed.
arma_simulate <- function(n, ar = numeric(0), ma = numeric(0)) {
  ss <- arma_state_space(ar, ma)
  ev <- eigen(ss$P0, symmetric = TRUE)
  a <- as.numeric(ev$vectors %*% (sqrt(pmax(ev$values, 0)) * stats::rnorm(ss$m)))
  e <- stats::rnorm(n)
  x <- numeric(n)
  x[1L] <- a[1L]
  for (t in seq_len(n - 1L) + 1L) {
    a <- ss$phi * a[1L] + c(a[-1L], 0) + ss$R * e[t]
    x[t] <- a[1L]
  }
  x
}

# Kalman filter for the zero-mean, unit-variance Gaussian ARMA process, in
# Harvey's state-space form. The innovation variance is 1 / gamma_0, so the
# stationary variance of x_t is exactly one: this is the constraint that makes
# x_t = qnorm(u_t) the normal scores of a copula. Missing values (NA) are
# skipped. Once the prediction covariance has converged the filter switches to
# the steady-state recursion (until the next missing value).
#
# Returns a list with the one-step predictive means 'mu' and variances 'F',
# the innovations 'v' (NA where x is missing), the Gaussian log-likelihood
# 'loglik' of the observed values, and the predictive state 'a' and covariance
# 'P' for time n + 1.
arma_filter <- function(x, ar = numeric(0), ma = numeric(0), tol = 1e-12) {
  n <- length(x)
  ss <- arma_state_space(ar, ma)
  m <- ss$m
  phi <- ss$phi
  Tm <- ss$Tm
  RR <- ss$RR
  Tstep <- function(a) phi * a[1L] + c(a[-1L], 0)

  a <- numeric(m)
  P <- ss$P0
  mu <- numeric(n)
  Fv <- numeric(n)
  v <- rep(NA_real_, n)
  ll <- 0
  steady <- FALSE
  K <- NULL
  for (t in seq_len(n)) {
    Ft <- P[1L, 1L]
    mu[t] <- a[1L]
    Fv[t] <- Ft
    if (is.na(x[t])) {
      a <- Tstep(a)
      P <- Tm %*% P %*% t(Tm) + RR
      steady <- FALSE
    } else {
      vt <- x[t] - a[1L]
      v[t] <- vt
      ll <- ll - 0.5 * (log(2 * pi) + log(Ft) + vt^2 / Ft)
      if (!steady) K <- P[, 1L] / Ft
      a <- Tstep(a + K * vt)
      if (!steady) {
        Pn <- Tm %*% (P - tcrossprod(P[, 1L]) / Ft) %*% t(Tm) + RR
        steady <- max(abs(Pn - P)) < tol
        P <- Pn
      }
    }
  }
  list(mu = mu, F = Fv, v = v, loglik = ll, a = a, P = P)
}

# Log-likelihood of copula data u in (0, 1) under the ARMA copula with the
# given (expanded) coefficients: the Gaussian log-likelihood of the normal
# scores minus the log standard normal density of the scores.
arma_copula_loglik <- function(u, ar = numeric(0), ma = numeric(0)) {
  x <- stats::qnorm(u)
  arma_filter(x, ar, ma)$loglik - sum(stats::dnorm(x, log = TRUE), na.rm = TRUE)
}
