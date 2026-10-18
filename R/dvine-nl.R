# Nonlinear (v2p) v-transforms in inverse-v-transformed (IVT) pair copulas: pure-R reference.
#
# See docs/sdvinecop-nonlinear-vt-spec.md. An "ivtnl" pair copula has the base copula `base` (as for
# "ivt"), the fulcrums `delta = c(delta1, delta2)` and the shapes `kappa = c(kappa1, kappa2)` of the
# transforms T_1 = v2p(delta1, kappa1) of the earlier and T_2 = v2p(delta2, kappa2) of the later
# variable of the pair. The v2p transform is UDP (T(U) is uniform when U is), so the density is
#   c_U(a, b) = c_V(T_1(a), T_2(b)),
# and the conditional distribution functions are integrals of that density over the variable that is
# being conditioned out:  P(U_1 <= a | U_2 = b) = int_0^a c_V(T_1(x), T_2(b)) dx,  and likewise for
# the later variable. The integral is split at the fulcrum and evaluated by graded Gauss-Legendre
# quadrature (24 nodes, grading 3, agreed 2026-10-18). When the transform of the integrated
# variable is linear (kappa = 1) the conditional distribution function is closed form (the base
# h-function at the other variable's transformed value), which is used exactly.

.dv_nl <- list(nodes = 24L, grade = 3)

dv_is_ivtnl <- function(pair) identical(pair$family, "ivtnl")

# v2p transform (the formulas of udp::v2p), vectorised; kappa = 1 is the linear transform.
dv_v2p <- function(x, d, k) {
  low <- 1 - x - (1 - d) * (x / d)^k
  upp <- x - d * ((1 - x) / (1 - d))^(1 / k)
  pmax(ifelse(x <= d, low, upp), 0)
}

# Derivative of the v2p transform (negative on the lower branch, positive on the upper one).
dv_v2p_slope <- function(x, d, k) {
  ifelse(x <= d, -1 - (1 - d) * k * x^(k - 1) / d^k, 1 + d * ((1 - x) / (1 - d))^(1 / k) / (k * (1 - x)))
}

# Resolved nonlinear IVT pair copula.
dv_ivtnl_pair <- function(base, tau, vt, kappa) {
  list(family = "ivtnl", rotation = 0, parameters = numeric(0), base = dv_ivt_base_pair(base, tau),
       delta = as.numeric(vt), kappa = as.numeric(kappa))
}

# Number of nonlinear lags: the last lag with |tau_k| above `nltau` (at least 1, at most K = length(tau)).
dv_nl_lags <- function(tau, nltau = 0.02) {
  K <- length(tau)
  above <- which(abs(tau) > nltau)
  min(max(1L, if (length(above)) max(above) else 0L), K)
}

# Pair copulas of the lags 1..K for the kpacf tau: lags up to L nonlinear, the rest linear IVT with the
# same fulcrums. `L = "auto"` applies the nltau rule.
dv_nl_pairs <- function(base, tau, vt, kappa, L = "auto", nltau = 0.02) {
  K <- length(tau)
  if (identical(L, "auto")) L <- dv_nl_lags(tau, nltau)
  L <- as.integer(min(max(L, 1L), K))
  lapply(seq_len(K), function(j) {
    if (j <= L) dv_ivtnl_pair(base, max(tau[j], 0), vt, kappa) else dv_ivt_pair(base, max(tau[j], 0), vt)
  })
}

# P(X <= x | Y = g) for a pair of IVT variables with transforms (dx, kx) of X and (dy, ky) of Y;
# x_first is TRUE when X is the earlier variable of the pair (h(a | b)) and FALSE when it is the
# later one (h(b | a)). `nodes` and `grade` set the quadrature.
dv_nl_cdf <- function(base, dx, kx, dy, ky, x, g, x_first, nodes = .dv_nl$nodes, grade = .dv_nl$grade) {
  n <- max(length(x), length(g))
  x <- rep_len(x, n)
  g <- rep_len(g, n)
  ty <- dv_clamp(dv_v2p(g, dy, ky))
  if (kx == 1) {
    # linear transform of the integrated variable: closed form
    tx <- dv_vlin(x, dx)
    hb <- if (x_first) dv_h(base, tx, ty, 2L) else dv_h(base, ty, tx, 1L)
    return(dx - dv_alin(x, dx) * hb)
  }
  gl <- gl_nodes(nodes)
  y <- gl$s
  w <- gl$w
  dens <- function(s, t2) {
    s <- dv_clamp(dv_v2p(as.vector(s), dx, kx))
    if (x_first) dv_dens(base, s, rep(t2, times = nodes)) else dv_dens(base, rep(t2, times = nodes), s)
  }
  m1 <- pmin(x, dx) # lower panel [0, min(x, dx)], nodes clustered at 0
  out <- rowSums(matrix(dens(outer(m1, y^grade), ty), n) * outer(m1, grade * y^(grade - 1) * w))
  up <- x > dx # upper panel [dx, x], nodes clustered at x
  if (any(up)) {
    len <- x[up] - dx
    f2 <- matrix(dens(dx + outer(len, 1 - (1 - y)^grade), ty[up]), sum(up))
    out[up] <- out[up] + rowSums(f2 * outer(len, grade * (1 - y)^(grade - 1) * w))
  }
  out
}

dv_ivtnl_dens <- function(pr, a, b) {
  dv_dens(pr$base, dv_clamp(dv_v2p(a, pr$delta[1L], pr$kappa[1L])), dv_clamp(dv_v2p(b, pr$delta[2L], pr$kappa[2L])))
}

dv_ivtnl_h <- function(pr, a, b, cond_var) {
  if (cond_var == 1L) dv_nl_cdf(pr$base, pr$delta[2L], pr$kappa[2L], pr$delta[1L], pr$kappa[1L], b, a, FALSE)
  else dv_nl_cdf(pr$base, pr$delta[1L], pr$kappa[1L], pr$delta[2L], pr$kappa[2L], a, b, TRUE)
}

# Inverse of h(b | a) in b: safeguarded Newton iteration on the monotone distribution function, with the
# density c_U(a, b) as derivative and bisection whenever a step leaves the bracket.
dv_ivtnl_hinv <- function(pr, a, p, tol = 1e-13, maxit = 80L) {
  n <- max(length(a), length(p))
  a <- rep_len(a, n)
  p <- rep_len(p, n)
  lo <- rep(0, n)
  hi <- rep(1, n)
  b <- dv_clamp(p)
  for (it in seq_len(maxit)) {
    f <- dv_ivtnl_h(pr, a, b, 1L) - p
    hi <- ifelse(f > 0, b, hi)
    lo <- ifelse(f <= 0, b, lo)
    if (all(abs(f) < tol | hi - lo < tol)) break
    d <- dv_ivtnl_dens(pr, a, b)
    bn <- b - f / d
    bad <- !is.finite(bn) | bn <= lo | bn >= hi
    b <- ifelse(bad, (lo + hi) / 2, bn)
  }
  b
}
