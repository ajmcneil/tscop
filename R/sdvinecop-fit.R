# Maximum likelihood fitting of tied sdvinecop models.

#' Fitted d-vine copula processes
#'
#' Class of objects returned by [tsfit()] for [sdvinecop-class] models.
#'
#' @slot model the fitted [sdvinecop-class] model (with its effective lag).
#' @slot data the copula data to which the model was fitted.
#' @slot pseudo logical: were the data pseudo-observations?
#' @slot estimate named numeric vector of estimates: the ARMA coefficients
#'   followed by the auxiliary parameters (`nu<lag>`, `delta<lag>`).
#' @slot se standard errors matching `estimate`; `NA` when not computed.
#' @slot vcov the matching covariance matrix; `NA` when not computed.
#' @slot se_method `"none"`, `"hessian"` or `"bootstrap"`.
#' @slot boot matrix of bootstrap estimates (no rows unless
#'   `se = "bootstrap"`).
#' @slot loglik maximised copula log-likelihood under the final effective
#'   lag.
#' @slot nobs number of observations.
#' @slot npar number of estimated parameters.
#' @slot convergence convergence code of the last optimisation (0 is
#'   success).
#' @slot message optimiser message.
#' @slot start the [sdvinecop-class] model used as starting values.
#' @slot Kstart the effective lag at the starting values.
#' @slot nrefit number of refits triggered by a change in the effective lag.
#' @slot dloglik the change in log-likelihood at the fitted parameters when
#'   the effective lag was updated at the end of the fit.
#' @slot template the original input if it was a `ts`, `zoo` or `xts` object,
#'   otherwise `NULL`.
#' @slot path for a model selected by [autodvinecop()], a data frame with one
#'   row per lag examined; an empty data frame otherwise.
#'
#' @exportClass sdvinecopfit
setClass("sdvinecopfit", slots = c(
  model = "sdvinecop", data = "numeric", pseudo = "logical",
  estimate = "numeric", se = "numeric", vcov = "matrix", se_method = "character",
  boot = "matrix", loglik = "numeric", nobs = "integer", npar = "integer",
  convergence = "integer", message = "character", start = "sdvinecop",
  Kstart = "integer", nrefit = "integer", dloglik = "numeric", template = "ANY",
  path = "data.frame"
))

# ---- parameter layout ----------------------------------------------------

# Free-scale layout: the unconstrained vector of the kpacf source (a Gaussian
# copula process, through its gc_free_ops closures) followed by one value per
# auxiliary parameter (t: nu = 2 + exp(z); BB1: delta = 1 + exp(z)).
dv_layout <- function(object) {
  ops <- gc_free_ops(object@kpacf)
  aux_lags <- which(vapply(object@pairs, function(sp) sp@family %in% .pc_aux_families, logical(1)))
  list(
    ops = ops, nk = ops$nfree,
    K0 = length(object@pairs), aux_lags = aux_lags,
    aux_fam = vapply(object@pairs[aux_lags], function(sp) sp@family, character(1)),
    explicit = lapply(object@pairs, function(sp) list(family = sp@family, pos = sp@pos, neg = sp@neg,
                                                      base = sp@base)),
    tail = object@tail,
    has_ivt = dv_has_ivt(object), vtpar = object@vtpar, vt_fixed = object@vt,
    nvt = if (!dv_has_ivt(object)) 0L else switch(object@vtpar, both = 2L, common = 1L, none = 0L)
  )
}

# Index of the v-transform parameters in the free vector (after the kpacf and
# auxiliary parameters), their values, and the free-scale image of fulcrums.
dv_vt_index <- function(lay) lay$nk + length(lay$aux_lags) + seq_len(lay$nvt)

dv_vt_from_free <- function(z, lay) {
  if (!lay$has_ivt) return(NULL)
  zv <- z[dv_vt_index(lay)]
  switch(lay$vtpar, both = stats::plogis(zv), common = rep(stats::plogis(zv), 2L), none = unname(lay$vt_fixed))
}

dv_vt_to_free <- function(vt, vtpar) {
  clamp <- function(x) pmin(pmax(x, 1e-3), 1 - 1e-3)
  switch(vtpar, both = stats::qlogis(clamp(vt)), common = stats::qlogis(clamp(vt[1L])), none = numeric(0))
}

dv_aux_from_free <- function(z, fam) ifelse(fam == "t", 2 + exp(z), 1 + exp(z))
dv_aux_to_free <- function(aux, fam) log(pmax(aux - ifelse(fam == "t", 2, 1), 1e-3))

# Natural coefficient vector (kpacf coefficients, then auxiliary parameters)
# from the free vector.
dv_free_to_coef <- function(z, lay) {
  na <- lay$nk
  c(lay$ops$coef(z[seq_len(na)]),
    dv_aux_from_free(z[na + seq_along(lay$aux_lags)], lay$aux_fam),
    if (lay$nvt > 0L) stats::plogis(z[dv_vt_index(lay)]))
}

# Resolved pair copulas of length K at the free vector z; NULL if infeasible.
dv_make_pairs <- function(z, lay, K) {
  na <- lay$nk
  tau <- (2 / pi) * asin(lay$ops$pacf(z[seq_len(na)], K))
  aux <- dv_aux_from_free(z[na + seq_along(lay$aux_lags)], lay$aux_fam)
  dv_make_pairs_tau(tau, aux, dv_vt_from_free(z, lay), lay, K)
}

# The same from the Kendall's tau values of lags 1..K (and the auxiliary
# parameters and fulcrums).
dv_make_pairs_tau <- function(tau, aux, vt, lay, K) {
  K0 <- lay$K0
  pairs <- vector("list", K)
  for (j in seq_len(min(K0, K))) {
    ex <- lay$explicit[[j]]
    a <- NA_real_
    ia <- match(j, lay$aux_lags)
    if (!is.na(ia)) a <- aux[ia]
    pr <- dv_resolve_vals(ex$family, ex$pos, ex$neg, a, tau[j], ex$base, vt)
    if (is.null(pr)) return(NULL)
    pairs[[j]] <- pr
  }
  if (K > K0) {
    tl <- dv_resolve_tail(lay$tail, tau[(K0 + 1L):K], vt)
    if (is.null(tl)) return(NULL)
    pairs[(K0 + 1L):K] <- tl
  }
  pairs
}

# Negative copula log-likelihood on the free scale at fixed effective lag K;
# a large finite value outside the feasible region or on numerical failure.
dv_negll <- function(z, lay, K, u) {
  na <- lay$nk
  if (any(!is.finite(z)) || (na > 0L && any(abs(z[seq_len(na)]) > lay$ops$zbound)) ||
      any(abs(z[na + seq_along(lay$aux_lags)]) > 12) ||
      (lay$nvt > 0L && any(abs(z[dv_vt_index(lay)]) > 6))) {
    return(1e10)
  }
  pairs <- tryCatch(dv_make_pairs(z, lay, K), error = function(e) NULL)
  if (is.null(pairs)) return(1e10)
  ll <- tryCatch(dv_loglik(u, pairs), error = function(e) -Inf)
  if (is.finite(ll)) -ll else 1e10
}

# Gradient of dv_negll with respect to the free vector z from the C++ kernel
# (forward-mode derivatives with respect to the pair copula parameters and the
# fulcrums, chained to z through the free-scale maps). The maps from z to the
# Kendall's tau values and from tau to the copula parameters are differentiated
# by central differences (cheap, deterministic, smooth), so the whole data
# recursion is differentiated exactly. NULL when the kernel does not apply:
# unsupported families, an infeasible point, more than 12 free parameters, or a
# rotation that changes within the difference step; the caller then uses
# finite differences of dv_negll.
dv_negll_grad <- function(z, lay, K, u) {
  if (!isTRUE(getOption("tscop.cpp", TRUE)) || !isTRUE(getOption("tscop.cppgrad", TRUE))) return(NULL)
  na <- lay$nk
  P <- length(z)
  if (P < 1L || P > 12L || length(lay$aux_lags) > 0L) return(NULL)
  if (any(!is.finite(z)) || (na > 0L && any(abs(z[seq_len(na)]) > lay$ops$zbound)) ||
      (lay$nvt > 0L && any(abs(z[dv_vt_index(lay)]) > 6))) return(NULL)
  zk <- z[seq_len(na)]
  tau_of <- function(zz) (2 / pi) * asin(lay$ops$pacf(zz, K))
  tau <- tau_of(zk)
  vt <- dv_vt_from_free(z, lay)
  pairs <- dv_make_pairs_tau(tau, numeric(0), vt, lay, K)
  if (is.null(pairs)) return(NULL)
  e <- dv_encode(pairs)
  if (is.null(e)) return(NULL)
  # d tau / d z (K x na)
  dtau <- matrix(0, K, P)
  for (k in seq_len(na)) {
    h <- 1e-6 * max(1, abs(zk[k]))
    zp <- zk; zp[k] <- zp[k] + h
    zm <- zk; zm[k] <- zm[k] - h
    dtau[, k] <- (tau_of(zp) - tau_of(zm)) / (2 * h)
  }
  # d parameter / d tau, lag by lag
  # (the numerical inverse of rvinecopulib for Joe and Frank has noise of about
  # 3e-10 in tau, so the step is 1e-4, giving a relative slope error of about 1e-5
  # at worst, and the truncation error of the central difference is far smaller)
  ht <- pmin(1e-4, 0.5 * abs(tau), 0.5 * (1 - abs(tau)))
  ht[ht < 1e-12] <- 1e-12
  pp <- dv_make_pairs_tau(tau + ht, numeric(0), vt, lay, K)
  pm <- dv_make_pairs_tau(tau - ht, numeric(0), vt, lay, K)
  if (is.null(pp) || is.null(pm)) return(NULL)
  ep <- dv_encode(pp)
  em <- dv_encode(pm)
  if (is.null(ep) || is.null(em)) return(NULL)
  # lags whose rotation or family changes within the step have |tau| < 2e-12: their
  # derivative is set to zero (their dependence on the parameters is negligible)
  same <- ep$rot == e$rot & em$rot == e$rot & ep$fam == e$fam & em$fam == e$fam &
    ep$brot == e$brot & em$brot == e$brot & ep$bfam == e$bfam & em$bfam == e$bfam
  dpar_dtau <- ifelse(same, (ep$par - em$par) / (2 * ht), 0)
  dbpar_dtau <- ifelse(same, (ep$bpar - em$bpar) / (2 * ht), 0)
  dpar <- dtau * dpar_dtau
  dbpar <- dtau * dbpar_dtau
  dd1 <- dd2 <- numeric(P)
  if (lay$nvt > 0L) {
    iv <- dv_vt_index(lay)
    if (lay$vtpar == "both") {
      dd1[iv[1L]] <- e$d1 * (1 - e$d1)
      dd2[iv[2L]] <- e$d2 * (1 - e$d2)
    } else {  # common
      dd1[iv[1L]] <- e$d1 * (1 - e$d1)
      dd2[iv[1L]] <- e$d2 * (1 - e$d2)
    }
  }
  g <- dv_loglik_grad_cpp(u, e$fam, e$rot, e$par, dpar, e$ivt, e$bfam, e$brot, e$bpar, dbpar,
                          e$d1, e$d2, dd1, dd2)
  if (!is.finite(g$value) || any(!is.finite(g$gradient))) return(NULL)
  -g$gradient
}

# Gradient for optim(): the exact gradient where available, otherwise central
# finite differences of dv_negll.
dv_negll_gr <- function(z, lay, K, u) {
  g <- dv_negll_grad(z, lay, K, u)
  if (!is.null(g)) return(g)
  vapply(seq_along(z), function(k) {
    h <- 1e-6
    zp <- z; zp[k] <- zp[k] + h
    zm <- z; zm[k] <- zm[k] - h
    (dv_negll(zp, lay, K, u) - dv_negll(zm, lay, K, u)) / (2 * h)
  }, numeric(1))
}

# Effective lag at the free vector z.
dv_K_at <- function(z, lay, tautol, maxlag) {
  lay$ops$effective_lag(z[seq_len(lay$nk)], lay$K0, tautol, maxlag)
}

# ---- effective-lag rule ----------------------------------------------------

# Wrap a fixed-K optimisation in the rule of the spec: fit with K held fixed;
# at the optimum recompute K and the change in log-likelihood; accept if the
# change is below 'tol' (storing the new K, or keeping the old one when the new
# is smaller and the change is not negligible); otherwise refit with the new K.
# K never decreases within a fit; at most 'max_refit' refits. fit_fixed(K, z)
# returns list(z, ...); K_of(z) the effective lag; ll_of(z, K) the
# log-likelihood.
dv_k_rule <- function(fit_fixed, K_of, ll_of, K_start, z0, tol = 0.01, max_refit = 3L) {
  K <- K_start
  z <- z0
  nrefit <- 0L
  repeat {
    fit <- fit_fixed(K, z)
    z <- fit$z
    Knew <- K_of(z)
    if (Knew == K) {
      return(list(fit = fit, z = z, K = K, Kstart = K_start, nrefit = nrefit, dloglik = 0, converged = TRUE))
    }
    d <- ll_of(z, Knew) - ll_of(z, K)
    if (Knew < K) {
      if (abs(d) < tol) {
        return(list(fit = fit, z = z, K = Knew, Kstart = K_start, nrefit = nrefit, dloglik = d, converged = TRUE))
      }
      return(list(fit = fit, z = z, K = K, Kstart = K_start, nrefit = nrefit, dloglik = 0, converged = TRUE))
    }
    if (abs(d) < tol) {
      return(list(fit = fit, z = z, K = Knew, Kstart = K_start, nrefit = nrefit, dloglik = d, converged = TRUE))
    }
    if (nrefit >= max_refit) {
      warning(sprintf("the effective lag was still changing after %d refits (now %d); the fit may be slightly inconsistent with it",
                      max_refit, Knew), call. = FALSE)
      return(list(fit = fit, z = z, K = Knew, Kstart = K_start, nrefit = nrefit, dloglik = d, converged = FALSE))
    }
    K <- Knew
    nrefit <- nrefit + 1L
  }
}

# ---- starting values -------------------------------------------------------

# Default starting values of the auxiliary parameters, given tau at their lags.
dv_default_aux <- function(spec_aux, fam, tau) {
  aux <- spec_aux
  for (i in seq_along(aux)) {
    if (!is.na(aux[i])) next
    if (fam[i] == "t") {
      aux[i] <- 10
    } else {
      dmax <- 1 / (1 - min(abs(tau[i]), 0.999))
      aux[i] <- min(1.5, 1 + 0.5 * (dmax - 1))
    }
  }
  aux
}

# Starting free vector, from init, from the model itself if it is fully
# specified, or from a Gaussian copula fit of the kpacf source.
dv_start <- function(object, init, u, control) {
  lay <- dv_layout(object)
  na <- lay$nk
  kp <- object@kpacf
  specaux <- vapply(object@pairs[lay$aux_lags], function(sp) sp@aux, numeric(1))
  vt <- if (lay$has_ivt) object@vt else numeric(0)
  if (!is.null(init)) {
    if (is(init, "sdvinecop")) {
      ok <- !is.null(init@kpacf) && identical(gc_structure(init@kpacf), gc_structure(kp)) &&
        length(init@pairs) == lay$K0 && identical(init@mode, "tied") &&
        identical(lapply(init@pairs, function(sp) sp@family), lapply(object@pairs, function(sp) sp@family))
      if (!ok) stop("'init' must be a tied sdvinecop with the same kpacf structure and pair copula families as 'object'")
      if (!gc_is_specified(init@kpacf)) stop("'init' has unspecified kpacf coefficients")
      kp <- init@kpacf
      specaux <- vapply(init@pairs[lay$aux_lags], function(sp) sp@aux, numeric(1))
      if (lay$has_ivt) vt <- init@vt
    } else if (is.numeric(init) && length(init) == na + length(lay$aux_lags) + lay$nvt) {
      kp <- gc_from_coef(kp, as.numeric(init[seq_len(na)]))
      specaux <- as.numeric(init[na + seq_along(lay$aux_lags)])
      if (lay$nvt > 0L) vt <- rep(as.numeric(init[na + length(lay$aux_lags) + seq_len(lay$nvt)]), length.out = 2L)
    } else {
      stop("'init' must be a sdvinecop model or a numeric vector of coefficients")
    }
  } else if (!gc_is_specified(kp)) {
    # Inverse-v-transformed pair copulas have cross-shaped dependence with no
    # rank correlation: the Gaussian kpacf is fitted to the v-transformed data
    # (one common fulcrum, the mean of the starting fulcrums), as in the vt-D-vine.
    up <- if (lay$has_ivt) pmin(pmax(dv_vlin(u, mean(vt)), 1e-6), 1 - 1e-6) else u
    kp <- gc_prefit(kp, up, control)
  }
  tau <- gc_tau(kp, max(lay$K0, 1L))
  aux <- dv_default_aux(specaux, lay$aux_fam, tau[lay$aux_lags])
  z0 <- c(lay$ops$to_free(kp), dv_aux_to_free(aux, lay$aux_fam),
          if (lay$has_ivt) dv_vt_to_free(vt, lay$vtpar))
  list(z0 = z0, lay = lay, kpacf = kp, aux = aux)
}

# The fitted model: kpacf coefficients and auxiliary parameters from z, with
# the stored effective lag K.
dv_build_model <- function(object, z, lay, K, limit) {
  kp <- lay$ops$build(z[seq_len(lay$nk)])
  aux <- dv_aux_from_free(z[lay$nk + seq_along(lay$aux_lags)], lay$aux_fam)
  pairs <- object@pairs
  for (i in seq_along(lay$aux_lags)) {
    sp <- pairs[[lay$aux_lags[i]]]
    pairs[[lay$aux_lags[i]]] <- pc(sp@family, sp@pos, sp@neg, aux[i])
  }
  vt <- dv_vt_from_free(z, lay)
  if (!is.null(vt)) names(vt) <- c("delta1", "delta2")
  m <- sdvinecop(kp, pairs = pairs, tail = object@tail, tautol = object@tautol,
                maxlag = object@maxlag, vt = vt, vtpar = object@vtpar)
  m@K <- as.integer(K)
  m@Klimit <- limit
  m
}

# ---- fitting --------------------------------------------------------------

# Full fitting procedure on complete data u: starting values, the K rule and
# the optimisation. Returns the pieces needed to build a sdvinecopfit.
dv_fit_all <- function(u, object, init, control) {
  stageK <- control$stageK
  control$stageK <- NULL
  st <- dv_start(object, init, u, control)
  lay <- st$lay
  z0 <- st$z0
  tautol <- object@tautol
  maxlag <- object@maxlag
  K_start <- dv_K_at(z0, lay, tautol, maxlag)$K
  if (dv_negll(z0, lay, K_start, u) >= 1e10) {
    stop("the starting values are infeasible for these data (a BB1 copula needs |tau| >= 1 - 1/delta at its lag)")
  }
  ctrl <- modifyList(list(maxit = 500, reltol = 1e-10), control)
  fit_fixed <- function(K, z, reltol = ctrl$reltol) {
    npar <- length(z)
    if (npar == 0L) return(list(z = z, value = dv_negll(z, lay, K, u), convergence = 0L, message = ""))
    c1 <- modifyList(ctrl, list(ndeps = rep(1e-6, npar), reltol = reltol))
    opt <- stats::optim(z, dv_negll, gr = dv_negll_gr, lay = lay, K = K, u = u, method = "BFGS", control = c1)
    if (opt$convergence != 0L || opt$value >= 1e10) {
      alt <- stats::optim(opt$par, dv_negll, lay = lay, K = K, u = u, method = "Nelder-Mead",
                          control = list(maxit = 5000, reltol = 1e-10))
      if (alt$value < opt$value) opt <- alt
    }
    list(z = opt$par, value = opt$value, convergence = opt$convergence,
         message = if (is.null(opt$message)) "" else opt$message)
  }
  K_of <- function(z) dv_K_at(z, lay, tautol, maxlag)$K
  ll_of <- function(z, K) -dv_negll(z, lay, K, u)
  # Staged fit: a long effective lag makes every evaluation expensive, so the
  # first optimisation uses the first K1 pair copulas only (loose tolerance);
  # its result is the starting point of the fit at the full lag.
  K1 <- if (is.null(stageK)) 10L else as.integer(stageK)
  user_start <- !is.null(init) || gc_is_specified(object@kpacf) # the user's values are the starting point
  if (!user_start && K1 >= 1L && K_start > 2L * K1 && length(z0) > 0L) {
    st1 <- fit_fixed(K1, z0, reltol = 1e-7)
    if (st1$convergence == 0L && st1$value < 1e10 && dv_negll(st1$z, lay, K_start, u) < 1e10) z0 <- st1$z
  }
  res <- dv_k_rule(fit_fixed, K_of, ll_of, K_of(z0), z0)
  limit <- dv_K_at(res$z, lay, tautol, maxlag)$limit
  if (res$K != K_of(res$z)) limit <- "tol"
  c(res, list(lay = lay, st = st, limit = limit, loglik = ll_of(res$z, res$K), u = u))
}

#' @describeIn tsfit Maximum likelihood fit of a tied [sdvinecop-class] model.
#'
#' @param bootstart for `se = "bootstrap"` on a [sdvinecop-class] model: `"fitted"`
#'   (default) starts each replicate's refit from the fitted values with the
#'   effective lag held fixed, which is much cheaper; `"user"` repeats the whole
#'   procedure from the original starting rules, including the effective-lag
#'   rule.
#'
#' @details For a [sdvinecop-class] model the ARMA coefficients and the
#'   auxiliary parameters (degrees of freedom of `t` copulas, `delta` of BB1
#'   copulas) are estimated jointly, by BFGS on an unconstrained scale with a
#'   Nelder-Mead fallback. Starting values come from `init`, otherwise from
#'   `object` if its ARMA coefficients are specified, otherwise from a Gaussian
#'   [armacop-class] fit to the data. The effective lag `K` is computed at the
#'   starting values and held fixed during the optimisation; at the optimum it
#'   is recomputed, and if the log-likelihood changes by 0.01 or more the model
#'   is refitted with the larger `K` (at most three times; `K` never decreases
#'   within a fit). The fitted object records the starting `K`, the number of
#'   refits and the final log-likelihood change.
#'
#'   When the starting `K` exceeds 20 and the starting values are not supplied
#'   by the user (`init` or a specified `object`) the optimisation is staged: a first fit
#'   uses only the first ten pair copulas (loose tolerance) and supplies the
#'   starting values of the fit at the full lag, which reduces the number of
#'   expensive evaluations. `control = list(stageK = 0)` switches this off, and
#'   `stageK` sets the number of pair copulas of the first stage.
#'
#'   Free-mode models cannot be fitted. The data must be complete.
#'   `se = "bootstrap"` is recommended: the Hessian standard errors are
#'   unreliable where a rotated copula switches rotation (Kendall's tau near
#'   zero), at the boundary of a BB1 copula's feasible region, and for
#'   pseudo-observations.
#'
#' @return For a [sdvinecop-class] model, an object of class
#'   [sdvinecopfit-class].
#' @export
setMethod("tsfit", signature(U = "ANY", object = "sdvinecop"),
  function(U, object, pseudo = TRUE, se = FALSE, B = 200, init = NULL,
           control = list(), bootstart = "fitted", ...) {
    if (object@mode != "tied") {
      stop("free-mode d-vine models cannot be fitted; specify a tied model (an armacop with pc() specifications)")
    }
    if (isTRUE(se)) se <- "bootstrap"
    if (!identical(se, FALSE) && !(is.character(se) && length(se) == 1L &&
                                   se %in% c("hessian", "bootstrap"))) {
      stop("'se' must be FALSE, \"hessian\" or \"bootstrap\"")
    }
    if (identical(se, "bootstrap") &&
        (!is.numeric(B) || length(B) != 1L || B < 10 || B != round(B))) {
      stop("'B' must be a whole number of at least 10")
    }
    bootstart <- match.arg(bootstart, c("fitted", "user"))
    u <- check_complete_copula_data(U, min_n = 20L)
    n <- length(u)
    fit <- dv_fit_all(u, object, init, control)
    lay <- fit$lay
    z <- fit$z
    K <- fit$K
    model <- dv_build_model(object, z, lay, K, fit$limit)
    est <- coef(model)
    cn <- names(est)
    npar <- length(est)
    V <- matrix(NA_real_, npar, npar, dimnames = list(cn, cn))
    boot <- matrix(NA_real_, 0L, npar, dimnames = list(NULL, cn))
    se_method <- "none"
    if (identical(se, "hessian") && npar > 0L) {
      H <- num_hessian(function(v) dv_negll(v, lay, K, u), z)
      Vz <- tryCatch(solve(H), error = function(e) NULL)
      if (is.null(Vz) || any(!is.finite(Vz)) || any(diag(Vz) <= 0) ||
          any(eigen(H, symmetric = TRUE, only.values = TRUE)$values <= 0)) {
        warning("the information matrix is not positive definite; standard errors not available")
      } else {
        J <- num_jacobian(function(v) dv_free_to_coef(v, lay), z)
        Vn <- J %*% Vz %*% t(J)
        V[] <- (Vn + t(Vn)) / 2
        se_method <- "hessian"
      }
    } else if (identical(se, "bootstrap") && npar > 0L) {
      pairs_hat <- dv_make_pairs(z, lay, K)
      Ub <- dv_simulate(pairs_hat, matrix(stats::runif(n * B), n, B))
      boot <- matrix(NA_real_, B, npar, dimnames = list(NULL, cn))
      for (b in seq_len(B)) {
        ub <- Ub[, b]
        if (pseudo) ub <- rank(ub, ties.method = "random") / (n + 1)
        est_b <- tryCatch(suppressWarnings({
          if (bootstart == "fitted") {
            ctrl <- modifyList(list(maxit = 500, reltol = 1e-10, ndeps = rep(1e-6, length(z))), control)
            opt <- stats::optim(z, dv_negll, gr = dv_negll_gr, lay = lay, K = K, u = ub, method = "BFGS", control = ctrl)
            if (opt$convergence != 0L || opt$value >= 1e10) stop("no convergence")
            dv_free_to_coef(opt$par, lay)
          } else {
            fb <- dv_fit_all(ub, object, init, control)
            dv_free_to_coef(fb$z, lay)
          }
        }), error = function(e) NULL)
        if (!is.null(est_b)) boot[b, ] <- est_b
      }
      good <- stats::complete.cases(boot)
      if (sum(good) >= 2L) {
        V[] <- stats::cov(boot[good, , drop = FALSE])
        se_method <- "bootstrap"
        if (sum(good) < B) {
          warning(sprintf("%d of %d bootstrap refits failed", B - sum(good), B))
        }
      } else {
        warning("too few successful bootstrap refits; standard errors not available")
      }
    }
    se_vec <- if (se_method == "none") rep(NA_real_, npar) else sqrt(diag(V))
    names(se_vec) <- cn
    start_model <- tryCatch(dv_build_model(object, fit$st$z0, lay, fit$Kstart, "tol"),
                            error = function(e) object)
    new("sdvinecopfit",
      model = model, data = u, pseudo = pseudo, estimate = est, se = se_vec, vcov = V,
      se_method = se_method, boot = boot, loglik = fit$loglik, nobs = n, npar = as.integer(npar),
      convergence = as.integer(fit$fit$convergence), message = fit$fit$message,
      start = start_model, Kstart = as.integer(fit$Kstart), nrefit = as.integer(fit$nrefit),
      dloglik = fit$dloglik,
      template = if (inherits(U, c("ts", "zoo", "xts"))) U else NULL,
      path = data.frame()
    )
  }
)

# ---- methods for fitted objects -------------------------------------------

#' @rdname tscop-methods
#' @export
setMethod("coef", "sdvinecopfit", function(object) object@estimate)

#' @rdname tscop-methods
#' @export
setMethod("logLik", "sdvinecopfit", function(object) {
  structure(object@loglik, df = object@npar, nobs = object@nobs, class = "logLik")
})

#' @rdname tscop-methods
#' @export
setMethod("vcov", "sdvinecopfit", function(object, ...) object@vcov)

#' @rdname tscop-methods
#' @export
setMethod("show", "sdvinecopfit", function(object) {
  cat("object class: sdvinecopfit\n")
  m <- object@model
  cat("kpacf: ", gc_describe(m@kpacf), "\n", sep = "")
  lab <- dv_explicit_labels(m)
  if (length(lab) > 0L) cat(length(lab), " explicit pair copulas: ", paste(lab, collapse = ", "), "\n", sep = "")
  cat("tail family: ", m@tail@family, "\n", sep = "")
  cat("effective lag K = ", m@K, " (tautol = ", m@tautol, ")",
      if (object@nrefit > 0L) paste0("; started at ", object@Kstart, ", ", object@nrefit, " refit(s)"),
      "\n", sep = "")
  cat("estimates:\n")
  if (object@se_method == "none") {
    print(object@estimate)
  } else {
    print(cbind(estimate = object@estimate, se = object@se))
    cat("standard errors: ", object@se_method, "\n", sep = "")
  }
  cat("loglik: ", format(object@loglik, digits = 7), "  npar: ", object@npar,
      "  nobs: ", object@nobs, "\n", sep = "")
  if (object@convergence != 0L) {
    cat("warning: optimiser did not converge (code ", object@convergence, ")\n", sep = "")
  }
  if (nrow(object@path) > 0L) {
    cat("selection path (see object@path):\n")
    print(object@path[, c("lag", "choice", "tau", "criterion", "improvement", "accepted")],
          row.names = FALSE, digits = 4)
  }
  invisible(object)
})

#' @rdname tssim
#' @export
setMethod("tssim", "sdvinecopfit", function(object, n, ...) {
  if (missing(n)) n <- length(object@data)
  tssim(object@model, n)
})

#' @rdname ktau
#' @export
setMethod("ktau", "sdvinecopfit", function(object, lagmax, ...) {
  if (missing(lagmax)) lagmax <- 20
  ktau(object@model, lagmax)
})

dv_fit_data <- function(object, data) {
  if (!missing(data)) return(data)
  if (is.null(object@template)) object@data else object@template
}

#' @rdname resid-method
#' @export
setMethod("resid", "sdvinecopfit", function(object, data, scale = "uniform", ...) {
  if (missing(data)) data <- dv_fit_data(object)
  dv_resid(object@model, data, scale)
})

#' @rdname resid-method
#' @export
setMethod("residuals", "sdvinecopfit", function(object, data, scale = "uniform", ...) {
  if (missing(data)) data <- dv_fit_data(object)
  dv_resid(object@model, data, scale)
})

#' @rdname predict-method
#' @export
setMethod("predict", "sdvinecopfit", function(object, data, h = 1, type = "df", x = NULL, ...) {
  if (missing(data)) data <- object@data
  dv_predict(object@model, data, h, type, x)
})
