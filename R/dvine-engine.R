# Engine for stationary d-vine copula processes: effective lag, resolution of
# pair copula specifications to concrete pair copulas, closed-form Gaussian
# pieces and the lag-by-lag likelihood recursion. All internal.
#
# A resolved pair copula is list(family, rotation, parameters) in the
# rvinecopulib convention applied to the pair (earlier variable, later
# variable). Gaussian pairs are recognised by family "gauss"/"gaussian".

# ---- effective lag -----------------------------------------------------

# Effective maximum lag of the truncated sequence of pair copulas: the last
# lag at which the ARMA Kendall partial autocorrelation exceeds 'tautol' in
# absolute value, capped at 'maxlag' and floored at max(K0, 1), where K0 is
# the number of explicit pair copulas. For a pure AR tail (no MA part) the
# PACF vanishes beyond the AR order, so the lag is the AR order exactly.
# 'period' (NA if none) sets the safety margin for seasonal bumps. Returns
# list(K, limit) with limit one of "ar", "tol", "maxlag", "cap", "floor".
dv_effective_lag <- function(ar, ma, K0 = 0L, tautol = 1e-4, maxlag = Inf,
                             period = NA, cap = 5000L) {
  floorK <- max(as.integer(K0), 1L)
  if (length(ma) == 0L) {
    K <- length(ar)
    limit <- "ar"
  } else {
    margin <- max(10L, 2L * (if (is.na(period)) 0L else as.integer(period)))
    Lmax <- cap
    if (is.finite(maxlag)) Lmax <- min(cap, max(50L, as.integer(maxlag) + margin))
    L <- min(50L, Lmax)
    repeat {
      tau <- (2 / pi) * asin(arma_pacf(ar, ma, L))
      ex <- which(abs(tau) > tautol)
      last <- if (length(ex)) max(ex) else 0L
      if (L - last >= margin) {
        K <- last
        limit <- "tol"
        break
      }
      if (L >= Lmax) {
        K <- last
        limit <- "tol"
        if (!is.finite(maxlag) || last < maxlag) {
          warning(sprintf(
            "the Kendall partial autocorrelations are still above 'tautol' at lag %d; truncating there (increase 'tautol' or set 'maxlag')",
            L), call. = FALSE)
          limit <- "cap"
        }
        break
      }
      L <- min(2L * L, Lmax)
    }
  }
  if (is.finite(maxlag) && K > maxlag) {
    K <- as.integer(maxlag)
    limit <- "maxlag"
  }
  if (K < floorK) {
    K <- floorK
    limit <- "floor"
  }
  list(K = as.integer(K), limit = limit)
}

# ---- resolution of pair copula specifications --------------------------

dv_is_gauss <- function(pair) pair$family %in% c("gauss", "gaussian")

# Copula parameter of a one-parameter family at Kendall's tau (>= 0 for the
# families without negative dependence). rvinecopulib::ktau_to_par is not
# vectorised and returns a matrix, so closed forms are used where they exist.
dv_ktau_to_par <- function(family, tau) {
  switch(family,
    clayton = pmax(2 * tau / (1 - tau), 1e-10),
    gumbel = 1 / (1 - tau),
    vapply(tau, function(t) as.numeric(rvinecopulib::ktau_to_par(family, t)), numeric(1))
  )
}

# Resolve a pc specification at Kendall's tau 'tau' to a pair copula, or NULL
# if the specification is infeasible at that tau (BB1 only).
dv_resolve_pair <- function(spec, tau) {
  fam <- spec@family
  if (fam %in% .pc_aux_families && is.na(spec@aux)) {
    stop("the auxiliary parameter of the '", fam, "' specification is missing")
  }
  if (fam == "gauss") {
    return(list(family = "gauss", rotation = 0, parameters = sin(pi * tau / 2)))
  }
  if (fam == "t") {
    return(list(family = "t", rotation = 0, parameters = c(sin(pi * tau / 2), spec@aux)))
  }
  if (fam == "frank") {
    return(list(family = "frank", rotation = 0,
                parameters = if (tau == 0) 1e-10 else dv_ktau_to_par("frank", tau)))
  }
  rot <- if (tau >= 0) spec@pos else spec@neg
  at <- abs(tau)
  if (fam == "bb1") {
    theta_c <- dv_ktau_to_par("clayton", at)
    th <- (theta_c + 2) / spec@aux - 2
    if (!(th > 1e-10)) return(NULL)
    return(list(family = "bb1", rotation = rot, parameters = c(th, spec@aux)))
  }
  list(family = fam, rotation = rot, parameters = dv_ktau_to_par(fam, at))
}

# Convert a bicop_dist object to a resolved pair copula.
dv_from_bicop <- function(b) {
  list(family = b$family, rotation = b$rotation, parameters = as.numeric(b$parameters))
}

# Kendall's tau of a resolved pair copula.
dv_pair_ktau <- function(pair) {
  rvinecopulib::par_to_ktau(pair$family, pair$rotation, pair$parameters)
}

# ---- likelihood recursion ----------------------------------------------

# Log-density of the d-vine copula process truncated at K = length(pairs):
# sum over trees j = 1..K of the log pair-copula densities at the forward and
# backward Rosenblatt variables. u must be a complete vector in (0, 1). Runs of
# Gaussian pair copulas are evaluated on the normal-score scale, where the
# h-functions are linear, so no pnorm/qnorm calls are needed between them.
# Trees beyond n - 1 do not enter the density.
dv_loglik <- function(u, pairs) {
  n <- length(u)
  K <- min(length(pairs), n - 1L)
  if (K < 1L) return(0)
  A <- u[-n] # earlier variable of each consecutive pair
  B <- u[-1L] # later variable
  zscale <- FALSE
  ll <- 0
  for (j in seq_len(K)) {
    pr <- pairs[[j]]
    if (dv_is_gauss(pr)) {
      if (!zscale) {
        A <- stats::qnorm(A)
        B <- stats::qnorm(B)
        zscale <- TRUE
      }
      rho <- pr$parameters[1L]
      r2 <- rho * rho
      ll <- ll + sum(-0.5 * log1p(-r2) - (r2 * (A * A + B * B) - 2 * rho * A * B) / (2 * (1 - r2)))
      if (j < K) {
        m <- length(A)
        s <- sqrt(1 - r2)
        Anew <- (A[-m] - rho * B[-m]) / s
        Bnew <- (B[-1L] - rho * A[-1L]) / s
        A <- Anew
        B <- Bnew
      }
    } else {
      if (zscale) {
        A <- stats::pnorm(A)
        B <- stats::pnorm(B)
        zscale <- FALSE
      }
      d <- rvinecopulib::dbicop(cbind(A, B), pr$family, pr$rotation, pr$parameters)
      ll <- ll + sum(log(d))
      if (j < K) {
        m <- length(A)
        v <- cbind(A, B)
        Anew <- rvinecopulib::hbicop(v[-m, , drop = FALSE], 2L, pr$family, pr$rotation, pr$parameters)
        Bnew <- rvinecopulib::hbicop(v[-1L, , drop = FALSE], 1L, pr$family, pr$rotation, pr$parameters)
        A <- as.numeric(Anew)
        B <- as.numeric(Bnew)
      }
    }
  }
  if (is.finite(ll)) ll else -Inf
}
