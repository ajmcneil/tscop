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
dv_is_ast <- function(pair) identical(pair$family, "ast")

# Rotations (rvinecopulib convention, applied to the pair (earlier, later)) of
# a copula given by closures: d(a, b) the density, h(a, b, cond_var) the
# h-function (cond_var 1: P(U2 <= b | U1 = a); 2: P(U1 <= a | U2 = b)) and
# hinv(a, p) the inverse of h(a, ., 1) in its second argument. Used for copulas
# that rvinecopulib does not implement (ast); validated against rvinecopulib's
# own rotated families.
dv_rotate <- function(base, rotation) {
  d <- base$d
  h <- base$h
  hinv <- base$hinv
  switch(as.character(rotation),
    "0" = base,
    "180" = list(
      d = function(a, b) d(1 - a, 1 - b),
      h = function(a, b, cv) 1 - h(1 - a, 1 - b, cv),
      hinv = function(a, p) 1 - hinv(1 - a, 1 - p)
    ),
    "90" = list(
      d = function(a, b) d(1 - a, b),
      h = function(a, b, cv) if (cv == 1L) h(1 - a, b, 1L) else 1 - h(1 - a, b, 2L),
      hinv = function(a, p) hinv(1 - a, p)
    ),
    "270" = list(
      d = function(a, b) d(a, 1 - b),
      h = function(a, b, cv) if (cv == 1L) 1 - h(a, 1 - b, 1L) else h(a, 1 - b, 2L),
      hinv = function(a, p) 1 - hinv(a, 1 - p)
    ),
    stop("rotation must be 0, 90, 180 or 270")
  )
}

# Closures of an ast pair copula with its rotation.
dv_ast_closures <- function(pr) {
  ac <- udp::astcopula(pr$parameters[1L])
  base <- list(
    d = function(a, b) udp::dastcopula(a, b, ac),
    h = function(a, b, cv) udp::hastcopula(a, b, ac, cond_var = cv),
    hinv = function(a, p) udp::hastcopula(a, p, ac, cond_var = 1L, inverse = TRUE)
  )
  dv_rotate(base, pr$rotation)
}

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

# Resolve a pair copula (family, rotations for positive and negative
# dependence, auxiliary parameter) at Kendall's tau 'tau' to a resolved pair
# copula, or NULL if infeasible at that tau (BB1 only).
dv_resolve_vals <- function(fam, pos, neg, aux, tau) {
  if (fam %in% .pc_aux_families && is.na(aux)) {
    stop("the auxiliary parameter of the '", fam, "' specification is missing")
  }
  if (fam == "gauss") {
    return(list(family = "gauss", rotation = 0, parameters = sin(pi * tau / 2)))
  }
  if (fam == "t") {
    return(list(family = "t", rotation = 0, parameters = c(sin(pi * tau / 2), aux)))
  }
  if (fam == "frank") {
    return(list(family = "frank", rotation = 0,
                parameters = if (tau == 0) 1e-10 else dv_ktau_to_par("frank", tau)))
  }
  rot <- if (tau >= 0) pos else neg
  at <- abs(tau)
  if (fam == "ast") {
    if (at < 1e-9) return(list(family = "indep", rotation = 0, parameters = numeric(0)))
    return(list(family = "ast", rotation = rot,
                parameters = udp::astcopula_nu(min(at, 0.999999))))
  }
  if (fam == "bb1") {
    theta_c <- dv_ktau_to_par("clayton", at)
    th <- (theta_c + 2) / aux - 2
    if (!(th > 1e-10)) return(NULL)
    return(list(family = "bb1", rotation = rot, parameters = c(th, aux)))
  }
  list(family = fam, rotation = rot, parameters = dv_ktau_to_par(fam, at))
}

# Resolve a pc specification at Kendall's tau 'tau' (see dv_resolve_vals).
dv_resolve_pair <- function(spec, tau, aux = spec@aux) {
  dv_resolve_vals(spec@family, spec@pos, spec@neg, aux, tau)
}

# Resolve the tail family at a vector of Kendall's tau values; the tail
# families have no auxiliary parameter and are always feasible.
dv_resolve_tail <- function(spec, tau) {
  fam <- spec@family
  if (fam == "gauss") {
    rho <- sin(pi * tau / 2)
    return(lapply(rho, function(r) list(family = "gauss", rotation = 0, parameters = r)))
  }
  if (fam == "frank") {
    par <- vapply(tau, function(t) if (t == 0) 1e-10 else dv_ktau_to_par("frank", t), numeric(1))
    return(lapply(par, function(p) list(family = "frank", rotation = 0, parameters = p)))
  }
  rot <- ifelse(tau >= 0, spec@pos, spec@neg)
  if (fam == "ast") {
    at <- abs(tau)
    nu <- rep(NA_real_, length(at))
    big <- at >= 1e-9
    if (any(big)) nu[big] <- udp::astcopula_nu(pmin(at[big], 0.999999))
    return(Map(function(n, r) {
      if (is.na(n)) list(family = "indep", rotation = 0, parameters = numeric(0))
      else list(family = "ast", rotation = r, parameters = n)
    }, nu, rot))
  }
  par <- dv_ktau_to_par(fam, abs(tau))
  Map(function(p, r) list(family = fam, rotation = r, parameters = p), par, rot)
}

# Convert a bicop_dist object to a resolved pair copula.
dv_from_bicop <- function(b) {
  list(family = b$family, rotation = b$rotation, parameters = as.numeric(b$parameters))
}

# Kendall's tau of a resolved pair copula.
dv_pair_ktau <- function(pair) {
  if (dv_is_ast(pair)) {
    tau0 <- udp::astcopula_tau(pair$parameters[1L])
    return(if (pair$rotation %in% c(90, 270)) -tau0 else tau0)
  }
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
    } else if (dv_is_ast(pr)) {
      if (zscale) {
        A <- stats::pnorm(A)
        B <- stats::pnorm(B)
        zscale <- FALSE
      }
      ll <- ll + sum(log(dv_dens(pr, A, B)))
      if (j < K) {
        m <- length(A)
        Anew <- dv_h(pr, A[-m], B[-m], 2L)
        Bnew <- dv_h(pr, A[-1L], B[-1L], 1L)
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

# ---- pair-level building blocks (uniform scale, vectorised) -------------
#
# For a pair (a, b) = (earlier variable, later variable):
#   dv_h(pr, a, b, 1) = h(b | a),   dv_h(pr, a, b, 2) = h(a | b),
#   dv_hinv(pr, a, p) = b such that h(b | a) = p.
# Gaussian pairs use closed forms. Results are clamped to the open unit
# interval so that they can be fed to later trees.

dv_clamp <- function(x, eps = 1e-12) pmin(pmax(x, eps), 1 - eps)

dv_h <- function(pr, a, b, cond_var) {
  if (dv_is_ast(pr)) {
    out <- dv_ast_closures(pr)$h(a, b, cond_var)
  } else if (dv_is_gauss(pr)) {
    rho <- pr$parameters[1L]
    s <- sqrt(1 - rho * rho)
    out <- if (cond_var == 1L) {
      stats::pnorm((stats::qnorm(b) - rho * stats::qnorm(a)) / s)
    } else {
      stats::pnorm((stats::qnorm(a) - rho * stats::qnorm(b)) / s)
    }
  } else {
    out <- as.numeric(rvinecopulib::hbicop(cbind(a, b), cond_var, pr$family,
                                           pr$rotation, pr$parameters))
  }
  dv_clamp(out)
}

dv_hinv <- function(pr, a, p) {
  if (dv_is_ast(pr)) {
    out <- dv_ast_closures(pr)$hinv(a, p)
  } else if (dv_is_gauss(pr)) {
    rho <- pr$parameters[1L]
    s <- sqrt(1 - rho * rho)
    out <- stats::pnorm(rho * stats::qnorm(a) + s * stats::qnorm(p))
  } else {
    out <- as.numeric(rvinecopulib::hbicop(cbind(a, p), 1L, pr$family, pr$rotation,
                                           pr$parameters, inverse = TRUE))
  }
  dv_clamp(out)
}

dv_dens <- function(pr, a, b) {
  if (dv_is_ast(pr)) {
    as.numeric(dv_ast_closures(pr)$d(a, b))
  } else if (dv_is_gauss(pr)) {
    rho <- pr$parameters[1L]
    x <- stats::qnorm(a)
    y <- stats::qnorm(b)
    r2 <- rho * rho
    exp(-0.5 * log1p(-r2) - (r2 * (x * x + y * y) - 2 * rho * x * y) / (2 * (1 - r2)))
  } else {
    as.numeric(rvinecopulib::dbicop(cbind(a, b), pr$family, pr$rotation, pr$parameters))
  }
}

# ---- forward Rosenblatt transform ---------------------------------------

# Forward Rosenblatt transform of a series: r_t = P(U_t <= u_t | previous
# min(t - 1, K) values), iid U(0, 1) under the model. Uses the same level
# recursion as dv_loglik; the first row of each level supplies the values at
# t = 2, ..., K + 1 that the density recursion discards.
dv_rosenblatt <- function(pairs, u) {
  n <- length(u)
  K <- min(length(pairs), n - 1L)
  r <- u
  if (K < 1L) return(r)
  A <- u[-n]
  Bv <- u[-1L]
  Bfinal <- numeric(0)
  for (j in seq_len(K)) {
    pr <- pairs[[j]]
    m <- length(A)
    r[j + 1L] <- dv_h(pr, A[1L], Bv[1L], 1L)
    if (j < K) {
      Anew <- dv_h(pr, A[-m], Bv[-m], 2L)
      Bnew <- dv_h(pr, A[-1L], Bv[-1L], 1L)
      A <- Anew
      Bv <- Bnew
    } else if (m > 1L) {
      Bfinal <- dv_h(pr, A[-1L], Bv[-1L], 1L)
    }
  }
  if (length(Bfinal) > 0L) r[(K + 2L):n] <- Bfinal
  r
}

# ---- simulation ----------------------------------------------------------

# Exact stationary simulation from independent U(0,1) innovations w (an n x
# nrep matrix, one column per replicate): u_t = R^{-1}_{min(t-1,K)}(w_t; past).
# Backward variables B_i(t - i), needed to invert the forward chain, are
# carried from step to step, so the cost is O(nK) pair-copula calls
# vectorised over replicates. dv_rosenblatt applied to a column of the result
# returns the corresponding column of w.
dv_simulate <- function(pairs, w) {
  n <- nrow(w)
  nrep <- ncol(w)
  K <- length(pairs)
  u <- matrix(0, n, nrep)
  Bprev <- matrix(0, 0L, nrep) # row i + 1 holds B_i(t - 1 - i)
  for (t in seq_len(n)) {
    k <- min(t - 1L, K)
    x <- w[t, ]
    Fm <- matrix(0, max(k, 1L), nrep) # row i holds F_{i-1}(t)
    if (k > 0L) {
      for (j in k:1L) {
        x <- dv_hinv(pairs[[j]], Bprev[j, ], x)
        Fm[j, ] <- x
      }
    }
    u[t, ] <- x
    nb <- min(k, K - 1L)
    Bnew <- matrix(0, nb + 1L, nrep)
    Bnew[1L, ] <- x
    if (nb >= 1L) {
      for (i in seq_len(nb)) {
        Bnew[i + 1L, ] <- dv_h(pairs[[i]], Bprev[i, ], Fm[i, ], 2L)
      }
    }
    Bprev <- Bnew
  }
  u
}

# ---- one-step predictive distribution ----------------------------------

# Backward variables B_i(n - i), i = 0, ..., k - 1, of the last k values of
# the data (k = min(K, n)): B_i(n - i) = P(U_{n-i} <= u_{n-i} | u_{n-i+1},
# ..., u_n). Level recursion on the window of the last k values; the last row
# of each level (dropped by the density recursion) supplies B_i(n - i).
dv_backward_tail <- function(pairs, u, k) {
  n <- length(u)
  w <- u[(n - k + 1L):n]
  out <- numeric(k)
  out[1L] <- w[k]
  if (k >= 2L) {
    A <- w[-k]
    Bv <- w[-1L]
    for (i in seq_len(k - 1L)) {
      pr <- pairs[[i]]
      m <- length(A)
      out[i + 1L] <- dv_h(pr, A[m], Bv[m], 2L)
      if (i < k - 1L) {
        Anew <- dv_h(pr, A[-m], Bv[-m], 2L)
        Bnew <- dv_h(pr, A[-1L], Bv[-1L], 1L)
        A <- Anew
        Bv <- Bnew
      }
    }
  }
  out
}

# Predictive distribution of the next value given the data: type "df"
# (distribution function), "qf" (quantile function) or "dens" (density),
# evaluated at x.
dv_predictive <- function(pairs, u, x, type) {
  n <- length(u)
  K <- length(pairs)
  k <- min(K, n)
  if (k < 1L) {
    return(switch(type, df = x, qf = x, dens = rep(1, length(x))))
  }
  Bt <- dv_backward_tail(pairs, u, k)
  switch(type,
    df = {
      f <- x
      for (j in seq_len(k)) f <- dv_h(pairs[[j]], rep(Bt[j], length(f)), f, 1L)
      f
    },
    qf = {
      f <- x
      for (j in k:1L) f <- dv_hinv(pairs[[j]], rep(Bt[j], length(f)), f)
      f
    },
    dens = {
      f <- x
      d <- rep(1, length(x))
      for (j in seq_len(k)) {
        a <- rep(Bt[j], length(f))
        d <- d * dv_dens(pairs[[j]], a, f)
        f <- dv_h(pairs[[j]], a, f, 1L)
      }
      d
    }
  )
}
