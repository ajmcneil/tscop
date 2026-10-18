# Engine for stationary d-vine copula processes: effective lag, resolution of
# pair copula specifications to concrete pair copulas, closed-form Gaussian
# pieces and the lag-by-lag likelihood recursion. All internal.
#
# A resolved pair copula is list(family, rotation, parameters) in the
# rvinecopulib convention applied to the pair (earlier variable, later
# variable). Gaussian pairs are recognised by family "gauss"/"gaussian".

# ---- effective lag -----------------------------------------------------

# Apply the cap 'maxlag' and the floor max(K0, 1) to an effective lag K.
dv_finish_lag <- function(K, limit, K0, maxlag) {
  floorK <- max(as.integer(K0), 1L)
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

# Effective lag from a function pacf_fn(L) returning the partial
# autocorrelations at lags 1..L: the last lag at which the Kendall partial
# autocorrelation exceeds 'tautol' in absolute value. The PACF is computed in
# doubling blocks of lags (50, 100, ...) and the search stops when the last
# exceedance lies at least 'margin' lags inside the computed range; a hard cap
# on the number of lags gives a warning if reached.
dv_effective_lag_blocks <- function(pacf_fn, K0 = 0L, tautol = 1e-4, maxlag = Inf,
                                    margin = 10L, cap = 5000L) {
  Lmax <- cap
  if (is.finite(maxlag)) Lmax <- min(cap, max(50L, as.integer(maxlag) + margin))
  L <- min(50L, Lmax)
  repeat {
    tau <- (2 / pi) * asin(pacf_fn(L))
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
  dv_finish_lag(K, limit, K0, maxlag)
}

# Effective maximum lag of the truncated sequence of pair copulas of an ARMA
# kpacf: the last lag at which the Kendall partial autocorrelation exceeds
# 'tautol' in absolute value, capped at 'maxlag' and floored at max(K0, 1),
# where K0 is the number of explicit pair copulas. For a pure AR kpacf (no MA
# part) the PACF vanishes beyond the AR order, so the lag is the AR order
# exactly. 'period' (NA if none) sets the safety margin for seasonal bumps.
# Returns list(K, limit) with limit one of "ar", "tol", "maxlag", "cap",
# "floor".
dv_effective_lag <- function(ar, ma, K0 = 0L, tautol = 1e-4, maxlag = Inf,
                             period = NA, cap = 5000L) {
  if (length(ma) == 0L) return(dv_finish_lag(length(ar), "ar", K0, maxlag))
  margin <- max(10L, 2L * (if (is.na(period)) 0L else as.integer(period)))
  dv_effective_lag_blocks(function(L) arma_pacf(ar, ma, L), K0, tautol, maxlag, margin, cap)
}

# ---- resolution of pair copula specifications --------------------------

dv_is_gauss <- function(pair) pair$family %in% c("gauss", "gaussian")
dv_is_ast <- function(pair) identical(pair$family, "ast")
dv_is_ivt <- function(pair) identical(pair$family, "ivt")

# Linear v-transform with fulcrum d and the signed scale factor a(x) (d for
# x <= d, d - 1 otherwise) of the IVT formulas (supplement to Dias, Han and
# McNeil, Proposition S3).
dv_vlin <- function(x, d) abs(x - d) / (d + (1 - 2 * d) * (x > d))
dv_alin <- function(x, d) d - (x > d)

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
# families without negative dependence). rvinecopulib::ktau_to_par is a scalar
# numerical inversion (noise of about 3e-10 in tau), so closed forms or the
# accurate vectorised inverses of R/tau-inverse.R are used where they exist.
dv_ktau_to_par <- function(family, tau) {
  switch(family,
    clayton = pmax(2 * tau / (1 - tau), 1e-10),
    gumbel = 1 / (1 - tau),
    joe = joe_par_fast(tau),
    frank = frank_par_fast(tau),
    vapply(tau, function(t) as.numeric(rvinecopulib::ktau_to_par(family, t)), numeric(1))
  )
}

# Resolve a pair copula (family, rotations for positive and negative
# dependence, auxiliary parameter) at Kendall's tau 'tau' to a resolved pair
# copula, or NULL if infeasible at that tau (BB1 only).
dv_resolve_vals <- function(fam, pos, neg, aux, tau, base = NA_character_, vt = NULL) {
  if (fam == "ivt") {
    if (is.null(vt)) stop("inverse-v-transformed pair copulas need the v-transform parameters")
    if (tau < -1e-10) return(NULL) # the base copula has positive dependence only
    return(dv_ivt_pair(base, max(tau, 0), vt))
  }
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
dv_resolve_pair <- function(spec, tau, aux = spec@aux, vt = NULL) {
  dv_resolve_vals(spec@family, spec@pos, spec@neg, aux, tau, spec@base, vt)
}

# Resolved base copula of an IVT pair copula at Kendall's tau (>= 0); the
# independence copula for a negligible tau.
dv_ivt_base_pair <- function(base, tau) {
  if (tau < 1e-9) return(list(family = "indep", rotation = 0, parameters = numeric(0)))
  switch(base,
    ast = list(family = "ast", rotation = 0, parameters = udp::astcopula_nu(min(tau, 0.999999))),
    joe = list(family = "joe", rotation = 0, parameters = dv_ktau_to_par("joe", tau)),
    gumbel = list(family = "gumbel", rotation = 0, parameters = dv_ktau_to_par("gumbel", tau)),
    clayton180 = list(family = "clayton", rotation = 180, parameters = dv_ktau_to_par("clayton", tau)),
    frank = list(family = "frank", rotation = 0, parameters = dv_ktau_to_par("frank", tau)),
    gauss = list(family = "gauss", rotation = 0, parameters = sin(pi * tau / 2)),
    stop("unknown base copula '", base, "'")
  )
}

# Resolved IVT pair copula: the base pair copula and the two deltas.
dv_ivt_pair <- function(base, tau, vt) {
  list(family = "ivt", rotation = 0, parameters = numeric(0),
       base = dv_ivt_base_pair(base, tau), delta = as.numeric(vt))
}

# Resolve the tail family at a vector of Kendall's tau values; the tail
# families have no auxiliary parameter and are always feasible.
dv_resolve_tail <- function(spec, tau, vt = NULL) {
  fam <- spec@family
  if (fam == "ivt") {
    if (is.null(vt)) stop("inverse-v-transformed pair copulas need the v-transform parameters")
    if (any(tau < -1e-10)) return(NULL)
    return(lapply(pmax(tau, 0), function(t) dv_ivt_pair(spec@base, t, vt)))
  }
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
  if (dv_is_ivt(pair) || dv_is_ivtnl(pair)) return(dv_pair_ktau(pair$base)) # the base copula's tau
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
dv_loglik_r <- function(u, pairs) {
  n <- length(u)
  K <- min(length(pairs), n - 1L)
  if (K < 1L) return(0)
  A <- u[-n] # earlier variable of each consecutive pair
  B <- u[-1L] # later variable
  zscale <- FALSE
  bs <- FALSE # is the recursion on the base scale of a run of IVT levels?
  zb <- FALSE # ... and are the base-scale values normal scores (Gaussian base)?
  ll <- 0
  for (j in seq_len(K)) {
    pr <- pairs[[j]]
    if (bs && !dv_is_ivt(pr)) {
      # leave the base scale: u-scale values from the base-scale ones
      if (zb) {
        x <- stats::pnorm(x)
        y <- stats::pnorm(y)
        zb <- FALSE
      }
      A <- dv_clamp(d1 - sA * x)
      B <- dv_clamp(d2 - sB * y)
      bs <- FALSE
    }
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
    } else if (dv_is_ivt(pr)) {
      # IVT pairs: by (S10) of the supplement of Dias, Han and McNeil,
      # V(h(a, b)) = h*(V(a), V(b)), so a run of IVT levels is evaluated on the
      # base scale throughout: the data are v-transformed once on entry, the
      # density at each level is the base density at (x, y), and the next
      # level's arguments are the base h-functions. The branch factors sA, sB of
      # the u-scale values (preserved by the h-functions) are carried so that
      # the u-scale values can be recovered on exit.
      if (!bs) {
        if (zscale) {
          A <- stats::pnorm(A)
          B <- stats::pnorm(B)
          zscale <- FALSE
        }
        d1 <- pr$delta[1L]
        d2 <- pr$delta[2L]
        x <- dv_vlin(A, d1)
        y <- dv_vlin(B, d2)
        sA <- dv_alin(A, d1)
        sB <- dv_alin(B, d2)
        bs <- TRUE
        zb <- FALSE
      }
      # a Gaussian base is evaluated on the normal-score scale, where its
      # h-functions are linear, so a run of Gaussian-base levels needs no
      # pnorm/qnorm calls between the levels
      gb <- dv_is_gauss(pr$base)
      if (gb && !zb) {
        x <- stats::qnorm(dv_clamp(x))
        y <- stats::qnorm(dv_clamp(y))
        zb <- TRUE
      } else if (!gb && zb) {
        x <- stats::pnorm(x)
        y <- stats::pnorm(y)
        zb <- FALSE
      }
      if (gb) {
        rho <- pr$base$parameters[1L]
        r2 <- rho * rho
        ll <- ll + sum(-0.5 * log1p(-r2) - (r2 * (x * x + y * y) - 2 * rho * x * y) / (2 * (1 - r2)))
        if (j < K) {
          m <- length(x)
          s <- sqrt(1 - r2)
          xn <- (x[-m] - rho * y[-m]) / s
          yn <- (y[-1L] - rho * x[-1L]) / s
        }
      } else {
        ll <- ll + sum(log(dv_dens(pr$base, x, y)))
        if (j < K) {
          m <- length(x)
          xn <- dv_h(pr$base, x[-m], y[-m], 2L)
          yn <- dv_h(pr$base, x[-1L], y[-1L], 1L)
        }
      }
      if (j < K) {
        sA <- sA[-m]
        sB <- sB[-1L]
        x <- xn
        y <- yn
      }
    } else if (dv_is_ivtnl(pr) || dv_is_ast(pr)) {
      # nonlinear IVT levels are evaluated on the u scale (V(h(a, b)) = h*(V(a), V(b)) holds for linear
      # v-transforms only), like the ast family
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

# ---- compiled likelihood ------------------------------------------------
#
# dv_loglik() evaluates the recursion in C++ (src/dvine_kernel.cpp) when every
# pair copula is supported by the kernel (independence, Gauss, Clayton, Gumbel,
# Joe, Frank with any rotation, and IVT pairs with one of these bases) and
# otherwise falls back to the R reference implementation dv_loglik_r().
# options(tscop.cpp = FALSE) forces the R implementation.

.cpp_family_codes <- c(indep = 0L, gauss = 1L, gaussian = 1L, clayton = 2L, gumbel = 3L, joe = 4L, frank = 5L, bb1 = 6L, ast = 7L, t = 8L)

# Integer/double vectors describing the resolved pair copulas for the kernel,
# or NULL if some pair copula is not supported. The absolute spherical t (ast) copula (code 7,
# rotation 0 only, plain or as an IVT base) is supported by the likelihood kernel alone, so it
# is encoded only when `ast = TRUE`.
dv_encode <- function(pairs, ast = FALSE, nl = FALSE) {
  K <- length(pairs)
  fam <- rot <- ivt <- bfam <- brot <- nlv <- integer(K)
  par <- bpar <- numeric(K)
  d1 <- d2 <- 0.5
  k1 <- k2 <- 1
  for (j in seq_len(K)) {
    pr <- pairs[[j]]
    if (identical(pr$family, "ivtnl")) {
      # nonlinear IVT levels: supported by the kernels for the Clayton, Gumbel, Joe, Frank and ast bases when the
      # caller handles them (nl = TRUE); the Gaussian base, and callers that do not (ktauhat), use the R implementation
      if (!nl) return(NULL)
      b <- pr$base
      code <- .cpp_family_codes[b$family]
      if (is.na(code) || code == 1L || (code == 7L && !ast)) return(NULL)
      nlv[j] <- 1L
      bfam[j] <- code
      brot[j] <- as.integer(b$rotation)
      bpar[j] <- if (length(b$parameters)) b$parameters[1L] else 0
      d1 <- pr$delta[1L]
      d2 <- pr$delta[2L]
      k1 <- pr$kappa[1L]
      k2 <- pr$kappa[2L]
    } else if (identical(pr$family, "ivt")) {
      b <- pr$base
      code <- .cpp_family_codes[b$family]
      if (is.na(code) || (code == 7L && !ast)) return(NULL)
      ivt[j] <- 1L
      bfam[j] <- code
      brot[j] <- as.integer(b$rotation)
      bpar[j] <- if (length(b$parameters)) b$parameters[1L] else 0
      d1 <- pr$delta[1L]
      d2 <- pr$delta[2L]
    } else {
      code <- .cpp_family_codes[pr$family]
      if (is.na(code) || (code == 7L && (!ast || pr$rotation != 0))) return(NULL)
      fam[j] <- code
      rot[j] <- as.integer(pr$rotation)
      par[j] <- if (length(pr$parameters)) pr$parameters[1L] else 0
      # the second parameter of a plain BB1 level travels in the (otherwise unused) bpar slot
      if (code == 6L || code == 8L) bpar[j] <- pr$parameters[2L]   # BB1: delta; t: the degrees of freedom
    }
  }
  list(fam = fam, rot = rot, par = par, ivt = ivt, bfam = bfam, brot = brot, bpar = bpar, d1 = d1, d2 = d2,
       nl = nlv, k1 = k1, k2 = k2)
}

dv_loglik <- function(u, pairs) {
  if (isTRUE(getOption("tscop.cpp", TRUE)) && length(pairs) > 0L) {
    e <- dv_encode(pairs, ast = TRUE, nl = TRUE)
    if (!is.null(e)) {
      return(dv_loglik_cpp(u, e$fam, e$rot, e$par, e$ivt, e$bfam, e$brot, e$bpar, e$d1, e$d2, e$nl, e$k1, e$k2))
    }
  }
  dv_loglik_r(u, pairs)
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
  if (dv_is_ivtnl(pr)) {
    out <- dv_ivtnl_h(pr, a, b, cond_var)
  } else if (dv_is_ivt(pr)) {
    d1 <- pr$delta[1L]
    d2 <- pr$delta[2L]
    x <- dv_vlin(a, d1)
    y <- dv_vlin(b, d2)
    out <- if (cond_var == 1L) d2 - dv_alin(b, d2) * dv_h(pr$base, x, y, 1L)
           else d1 - dv_alin(a, d1) * dv_h(pr$base, x, y, 2L)
  } else if (dv_is_ast(pr)) {
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
  if (dv_is_ivtnl(pr)) {
    out <- dv_ivtnl_hinv(pr, a, p)
  } else if (dv_is_ivt(pr)) {
    d1 <- pr$delta[1L]
    d2 <- pr$delta[2L]
    out <- d2 - dv_alin(p, d2) * dv_hinv(pr$base, dv_vlin(a, d1), dv_vlin(p, d2))
  } else if (dv_is_ast(pr)) {
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
  if (dv_is_ivtnl(pr)) {
    dv_ivtnl_dens(pr, a, b)
  } else if (dv_is_ivt(pr)) {
    dv_dens(pr$base, dv_vlin(a, pr$delta[1L]), dv_vlin(b, pr$delta[2L]))
  } else if (dv_is_ast(pr)) {
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
dv_rosenblatt_r <- function(pairs, u) {
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
dv_simulate_r <- function(pairs, w) {
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

# Compiled versions where every pair copula is supported by the kernels (see
# dv_encode); the R implementations above are the reference.
dv_simulate <- function(pairs, w) {
  if (isTRUE(getOption("tscop.cpp", TRUE)) && length(pairs) > 0L) {
    e <- dv_encode(pairs, ast = TRUE, nl = TRUE)
    if (!is.null(e)) return(dv_simulate_cpp(w, e$fam, e$rot, e$par, e$ivt, e$bfam, e$brot, e$bpar, e$d1, e$d2, e$nl, e$k1, e$k2))
  }
  dv_simulate_r(pairs, w)
}

dv_rosenblatt <- function(pairs, u) {
  if (isTRUE(getOption("tscop.cpp", TRUE)) && length(pairs) > 0L) {
    e <- dv_encode(pairs, ast = TRUE, nl = TRUE)
    if (!is.null(e)) return(dv_rosenblatt_cpp(u, e$fam, e$rot, e$par, e$ivt, e$bfam, e$brot, e$bpar, e$d1, e$d2, e$nl, e$k1, e$k2))
  }
  dv_rosenblatt_r(pairs, u)
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

# Predictive distribution of the next value given the data (R reference): type "df"
# (distribution function), "qf" (quantile function) or "dens" (density),
# evaluated at x.
dv_predictive_r <- function(pairs, u, x, type) {
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

# Compiled one-step predictive distribution where every pair copula is supported
# by the kernels (dv_encode), otherwise the R reference above.
dv_predictive <- function(pairs, u, x, type) {
  if (isTRUE(getOption("tscop.cpp", TRUE)) && length(pairs) > 0L && length(u) > 0L) {
    e <- dv_encode(pairs, ast = TRUE, nl = TRUE)
    if (!is.null(e)) {
      code <- match(type, c("df", "qf", "dens"))
      return(dv_predictive_cpp(as.numeric(u), as.numeric(x), code, e$fam, e$rot, e$par, e$ivt, e$bfam, e$brot, e$bpar, e$d1, e$d2,
                               e$nl, e$k1, e$k2))
    }
  }
  dv_predictive_r(pairs, u, x, type)
}
