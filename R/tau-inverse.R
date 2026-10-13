# Closed-form Kendall's tau of the Joe and Frank copulas, their derivatives in
# the copula parameter, and vectorised inverses (parameter from tau), accurate
# to about 1e-14. rvinecopulib::ktau_to_par is a scalar numerical inversion
# with noise of about 3e-10 in tau and about 100 microseconds per call.

# ---- Joe ---------------------------------------------------------------------
#
# tau(theta) = 1 + (2 / (2 - theta)) (psi(2) - psi(2 / theta + 1)). With
# r = 1 - 2 / theta (in [-1, 1) for theta >= 1) the Taylor series of the
# digamma function at 2, psi(2 + w) - psi(2) = sum_n (-1)^(n+1) (zeta(n+1) - 1) w^n
# with w = -(theta - 2) / theta, gives the convergent, cancellation-free form
#   tau(theta) = 1 - (2 / theta) sum_{n >= 1} c_n r^(n-1),   c_n = zeta(n+1) - 1,
# with terms of size about (r / 2)^n.

# zeta(k) - 1 for k = 2, ..., 70 by Euler-Maclaurin summation
.zeta_minus_1 <- function(k, M = 25L) {
  B <- c(1 / 6, -1 / 30, 1 / 42, -1 / 30, 5 / 66, -691 / 2730, 7 / 6)  # B2, B4, ..., B14
  vapply(k, function(s) {
    m <- seq_len(M - 1L)
    val <- sum(m^(-s)) + M^(1 - s) / (s - 1) + 0.5 * M^(-s)
    # + sum_j B_2j / (2j)! * s (s+1) ... (s+2j-2) M^(-s-2j+1)
    term <- s * M^(-s - 1)
    fact <- 2
    for (j in seq_along(B)) {
      val <- val + B[j] / fact * term
      # next term: multiply the rising factorial by (s+2j-1)(s+2j) and M^-2
      term <- term * (s + 2 * j - 1) * (s + 2 * j) * M^(-2)
      fact <- fact * (2 * j + 1) * (2 * j + 2)
    }
    val - 1
  }, numeric(1))
}
.joe_c <- .zeta_minus_1(2:71)    # c_n, n = 1..70

joe_tau <- function(theta) {
  r <- 1 - 2 / theta
  n <- length(.joe_c)
  s <- rep(.joe_c[n], length(theta))
  for (i in (n - 1L):1L) s <- s * r + .joe_c[i]     # Horner: sum c_n r^(n-1)
  1 - 2 / theta * s
}

# d tau / d theta
joe_dtau <- function(theta) {
  r <- 1 - 2 / theta
  n <- length(.joe_c)
  s <- rep(.joe_c[n], length(theta))
  ds <- rep(0, length(theta))
  for (i in (n - 1L):1L) {
    ds <- ds * r + s              # derivative of the Horner polynomial in r
    s <- s * r + .joe_c[i]
  }
  # tau = 1 - (2/theta) S(r), r = 1 - 2/theta, dr/dtheta = 2/theta^2
  2 / theta^2 * s - 2 / theta * ds * 2 / theta^2
}

# ---- Frank -------------------------------------------------------------------
#
# tau(theta) = 1 - (4 / theta) (1 - D1(theta)), D1(theta) = (1/theta) int_0^theta t / (e^t - 1) dt.
# For |theta| <= 2 the Bernoulli series of 1 - D1 is used, beyond that
# D1 = (pi^2/6 - sum_m e^(-m theta) (theta / m + 1 / m^2)) / theta.

.bern <- c(1 / 6, -1 / 30, 1 / 42, -1 / 30, 5 / 66, -691 / 2730, 7 / 6, -3617 / 510, 43867 / 798,
           -174611 / 330, 854513 / 138, -236364091 / 2730, 8553103 / 6, -23749461029 / 870,
           8615841276005 / 14322)
# 1 - D1(theta) = theta / 4 - sum_k B_2k theta^(2k) / ((2k + 1) (2k)!)
.frank_coef <- local({
  k <- seq_along(.bern)
  c2k <- -.bern / ((2 * k + 1) * factorial(2 * k))     # coefficient of theta^(2k) in 1 - D1
  c2k
})

frank_tau <- function(theta) {
  a <- abs(theta)
  out <- numeric(length(a))
  small <- a <= 2
  if (any(small)) {
    t <- a[small]
    t2 <- t * t
    # tau = -4 sum_k c_2k theta^(2k-1)
    s <- rep(.frank_coef[length(.frank_coef)], length(t))
    for (i in (length(.frank_coef) - 1L):1L) s <- s * t2 + .frank_coef[i]
    out[small] <- -4 * s * t
  }
  if (any(!small)) {
    t <- a[!small]
    m <- seq_len(60L)
    tail <- vapply(t, function(x) sum(exp(-m * x) * (x / m + 1 / m^2)), numeric(1))
    D1 <- (pi^2 / 6 - tail) / t
    out[!small] <- 1 - 4 / t * (1 - D1)
  }
  sign(theta) * out
}

# d tau / d theta (even in theta)
frank_dtau <- function(theta) {
  a <- abs(theta)
  out <- numeric(length(a))
  small <- a <= 2
  if (any(small)) {
    t <- a[small]
    t2 <- t * t
    k <- seq_along(.frank_coef)
    # tau = -4 sum_k c_2k theta^(2k-1); derivative -4 sum_k (2k-1) c_2k theta^(2k-2)
    s <- rep((2 * k[length(k)] - 1) * .frank_coef[length(k)], length(t))
    for (i in (length(k) - 1L):1L) s <- s * t2 + (2 * i - 1) * .frank_coef[i]
    out[small] <- -4 * s
  }
  if (any(!small)) {
    t <- a[!small]
    m <- seq_len(60L)
    tail <- vapply(t, function(x) sum(exp(-m * x) * (x / m + 1 / m^2)), numeric(1))
    D1 <- (pi^2 / 6 - tail) / t
    out[!small] <- 4 / t^2 * (1 - 2 * D1) + 4 / (t * expm1(t))
  }
  out
}

# ---- inversion ---------------------------------------------------------------

# Safeguarded Newton for an increasing function f(x) = target on x > lo with
# derivative df; vectorised over targets.
.invert_increasing <- function(f, df, target, lo, hi, x0, tol = 1e-15, maxit = 100L) {   # tolerance relative to x - lo
  x <- x0
  lo0 <- lo
  for (it in seq_len(maxit)) {
    fx <- f(x) - target
    above <- fx > 0
    hi <- ifelse(above, x, hi)
    lo <- ifelse(above, lo, x)
    xn <- x - fx / df(x)
    bad <- !is.finite(xn) | xn <= lo | xn >= hi
    xn[bad] <- 0.5 * (lo[bad] + hi[bad])
    done <- abs(xn - x) <= tol * abs(xn - lo0)
    x <- xn
    if (all(done)) break
  }
  x
}

# Joe parameter at Kendall's tau in [0, 1); 1 for tau = 0.
joe_par <- function(tau) {
  out <- rep(1, length(tau))
  ok <- tau > 1e-14
  if (any(ok)) {
    tt <- tau[ok]
    hi <- pmax(3, 4 / (1 - tt))
    while (any(joe_tau(hi) < tt)) hi <- ifelse(joe_tau(hi) < tt, 2 * hi, hi)
    lo <- rep(1, length(tt))
    out[ok] <- .invert_increasing(joe_tau, joe_dtau, tt, lo, hi, x0 = 1 + 1.5 * tt / (1 - tt^2 * 0.5))
  }
  out
}

# Frank parameter at Kendall's tau in (-1, 1); +-1e-10 for |tau| below 1.1e-11.
frank_par <- function(tau) {
  a <- abs(tau)
  out <- ifelse(tau < 0, -1e-10, 1e-10)   # |theta| >= 1e-10 (tau of 1.1e-11 or less is independence)
  ok <- a >= 1.1e-11
  if (any(ok)) {
    tt <- a[ok]
    hi <- pmax(10, 8 / (1 - tt))
    while (any(frank_tau(hi) < tt)) hi <- ifelse(frank_tau(hi) < tt, 2 * hi, hi)
    lo <- rep(0, length(tt))
    out[ok] <- sign(tau[ok]) * .invert_increasing(frank_tau, frank_dtau, tt, lo, hi, x0 = pmin(9 * tt, 0.5 * hi))
  }
  out
}

# ---- compiled versions (src/tau_inverse.cpp), used by the model code ---------

joe_par_fast <- function(tau) joe_par_cpp(as.numeric(tau), .joe_c)
frank_par_fast <- function(tau) frank_par_cpp(as.numeric(tau), .frank_coef)
joe_tau_d <- function(theta) joe_tau_cpp(as.numeric(theta), .joe_c)
frank_tau_d <- function(theta) frank_tau_cpp(as.numeric(theta), .frank_coef)
