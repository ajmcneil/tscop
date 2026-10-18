#' Stationary d-vine copula processes
#'
#' Stationary d-vine copula process of finite order \eqn{K}, built from a
#' sequence of pair copulas \eqn{C_1, \ldots, C_K} (independence beyond lag
#' \eqn{K}).
#'
#' There are two modes.
#'
#' **Tied mode** (the only mode that can be fitted). A known Gaussian
#' ARMA/SARMA [armacop-class] model (or an [arfactcop-class] model) supplies the Kendall partial
#' autocorrelation function (kpacf) \eqn{\tau_k = (2/\pi)\arcsin\alpha_k},
#' with \eqn{\alpha_k} the partial autocorrelation of the ARMA process at lag
#' \eqn{k}. The pair copula at lag \eqn{k} is a member of a family chosen for
#' that lag, calibrated to have Kendall's tau \eqn{\tau_k}. Lags
#' \eqn{1, \ldots, K_0} have explicit [pc-class] specifications (any
#' families, including Gaussian); lags \eqn{K_0 + 1, \ldots, K} use the
#' `tail` family.
#'
#' **Free mode** (known models only). An explicit list of
#' [rvinecopulib::bicop_dist()] objects with their own parameters, optionally
#' followed by a Gaussian ARMA tail. The rotation of a `bicop_dist` is
#' interpreted by [rvinecopulib::bicop_dist()] applied to the pair (earlier
#' variable, later variable).
#'
#' The effective lag \eqn{K} is computed once when the object is created and
#' stored: it is the last lag at which \eqn{|\tau_k|} exceeds `tautol`, capped
#' at `maxlag`, and at least the number of explicit pair copulas; for a pure
#' AR tail it is the AR order. The fitted and simulated process is therefore
#' always a finite-order Markov chain.
#'
#' @slot name short description.
#' @slot mode `"tied"` or `"free"`.
#' @slot kpacf a Gaussian copula process ([gausscop-class]), or `NULL`.
#' @slot pairs list of explicit pair copulas: [pc-class] objects (tied mode)
#'   or `bicop_dist` objects (free mode).
#' @slot tail the tail family: a [pc-class] object (tied mode), or `NULL`.
#' @slot tautol,maxlag the tolerance and cap defining the effective lag.
#' @slot vt the fulcrums `c(delta1, delta2)` of the linear v-transforms
#'   of the inverse-v-transformed (`"ivt"`) pair copulas, shared by all of them
#'   (`delta1` applies to the earlier and `delta2` to the later variable of
#'   a pair); empty if the model has no such pair copula.
#' @slot vtpar which fulcrums are parameters: `"both"`, `"common"` (one
#'   parameter, equal fulcrums) or `"none"` (fixed).
#' @slot K the effective lag (`NA` for a model structure whose ARMA
#'   coefficients are not yet specified).
#' @slot Klimit what determined `K`: `"ar"` (pure AR order), `"tol"`,
#'   `"maxlag"`, `"floor"` (number of explicit pair copulas), `"cap"`
#'   (internal limit on the lags examined) or `"none"` (free mode without an
#'   ARMA tail).
#'
#' @references Bladt, M. and McNeil, A. J. (2022). Time series copula models using d-vines
#'   and v-transforms. *Econometrics and Statistics*, 24, 27-48.
#'   \doi{10.1016/j.ecosta.2021.07.004}
#'
#'   Bladt, M. and McNeil, A. J. (2022). Time series with infinite-order partial
#'   copula dependence. *Dependence Modeling*, 10(1), 87-107.
#'   \doi{10.1515/demo-2022-0105}
#'
#'   Bladt, M., Dias, A., Han, J. and McNeil, A. J. (2026). Semiparametric and
#'   parametric distributional forecasting of univariate time series using
#'   non-Gaussian ARMA models based on D-vines. *Canadian Journal of Statistics*,
#'   54(1), e70033. \doi{10.1002/cjs.70033}
#'
#'   Dias, A., Han, J. and McNeil, A. J. (2027). GARCH copulas, v-transforms and
#'   D-vines for stochastic volatility. *Journal of Multivariate Analysis*, 217,
#'   105695. \doi{10.1016/j.jmva.2026.105695}
#'
#' @exportClass sdvinecop
setClass("sdvinecop", slots = c(
  name = "character", mode = "character", kpacf = "gausscopOrNULL",
  pairs = "list", tail = "pcOrNULL", tautol = "numeric", maxlag = "numeric",
  K = "integer", Klimit = "character", vt = "numeric", vtpar = "character"
))

# Does the model contain an inverse-v-transformed pair copula?
dv_has_ivt <- function(object) {
  any(vapply(object@pairs, function(sp) is(sp, "pc") && sp@family == "ivt", logical(1))) ||
    (!is.null(object@tail) && object@tail@family == "ivt")
}

setValidity("sdvinecop", function(object) {
  if (!(object@mode %in% c("tied", "free"))) return("'mode' must be \"tied\" or \"free\"")
  if (length(object@tautol) != 1L || !(object@tautol > 0)) return("'tautol' must be a positive number")
  if (length(object@maxlag) != 1L || is.na(object@maxlag) || object@maxlag < 1) {
    return("'maxlag' must be at least 1 (use Inf for no cap)")
  }
  K0 <- length(object@pairs)
  if (object@mode == "tied") {
    if (is.null(object@kpacf)) return("a tied model needs a 'kpacf' component")
    if (!all(vapply(object@pairs, is, logical(1), "pc"))) return("tied 'pairs' must be pc() specifications")
    if (is.null(object@tail)) return("a tied model needs a tail family")
    if (!(object@tail@family %in% .pc_tail_families)) {
      return(paste0("the tail family must be one of ", paste(.pc_tail_families, collapse = ", "),
                    " (families that tend to independence as Kendall's tau tends to 0)"))
    }
  } else {
    if (!all(vapply(object@pairs, is, logical(1), "bicop_dist"))) {
      return("free 'pairs' must be bicop_dist objects")
    }
    if (K0 == 0L && is.null(object@kpacf)) return("specify 'pairs' or 'kpacf'")
  }
  if (!is.na(object@K) && object@K < max(K0, 1L)) return("'K' must be at least the number of explicit pairs")
  if (length(object@vtpar) != 1L || !(object@vtpar %in% c("both", "common", "none"))) {
    return("'vtpar' must be \"both\", \"common\" or \"none\"")
  }
  if (dv_has_ivt(object)) {
    if (length(object@vt) != 2L || anyNA(object@vt) || any(object@vt <= 0 | object@vt >= 1)) {
      return("the v-transform deltas 'vt' must be two numbers in (0, 1)")
    }
    if (object@vtpar == "common" && object@vt[1L] != object@vt[2L]) {
      return("with vtpar = \"common\" the two deltas must be equal")
    }
  } else if (length(object@vt) > 0L) {
    return("'vt' is only used with inverse-v-transformed pair copulas")
  }
  TRUE
})

#' Constructor for stationary d-vine copula processes
#'
#' Creates a [sdvinecop-class] object in tied mode (`pairs` are
#' [pc()] specifications) or free mode (`pairs` are
#' [rvinecopulib::bicop_dist()] objects).
#'
#' @param kpacf a Gaussian copula process ([armacop-class] or
#'   [arfactcop-class]). In tied mode its Kendall partial
#'   autocorrelations calibrate every pair copula; it may be a model structure
#'   (`armacop(order = ...)`) if the model is only to be passed to a fitting
#'   function. In free mode it is optional and supplies a Gaussian tail.
#' @param pairs list of explicit pair copulas for lags `1, ..., length(pairs)`
#'   (a single [pc()] specification is the list of length one; `rep(pc(...), 3)`
#'   repeats one specification over three lags):
#'   all [pc()] specifications (tied mode) or all `bicop_dist` objects (free
#'   mode). A `pc("gauss")` entry keeps a lag Gaussian.
#' @param tail the family used for all lags beyond `pairs` in tied mode:
#'   `"gauss"` (default), `"frank"`, `"clayton"`, `"gumbel"`, `"joe"` or
#'   `"ast"`, or a [pc()] specification of one of these, or of an `"ivt"`
#'   pair copula (`pc("ivt", base = ...)`, with base `"ast"`, `"joe"`,
#'   `"gumbel"`, `"clayton180"`, `"frank"` or `"gauss"`; it uses the fulcrums `vt` of the
#'   model; the `"gumbel"` and `"gauss"` bases make the fulcrum estimates
#'   unreliable, see [pc()]). Families whose Kendall's tau cannot tend to zero with their
#'   auxiliary parameters fixed (`"t"`, `"bb1"`) are rejected. Free mode has a
#'   Gaussian tail.
#' @param vt the fulcrums of the linear v-transforms of the `"ivt"` pair
#'   copulas, two numbers in (0, 1) (default `c(0.5, 0.5)` when the model has an
#'   `"ivt"` pair copula; an error otherwise). They are starting values when
#'   the model is fitted, unless `vtpar = "none"`.
#' @param vtpar `"both"` (default): both fulcrums are estimated; `"common"`:
#'   one fulcrum for both (the vt-D-vine of Bladt and McNeil, with `vt` equal in
#'   its two components); `"none"`: the fulcrums are fixed at `vt`.
#' @param tautol tolerance on the absolute Kendall's tau defining the
#'   effective lag.
#' @param maxlag cap on the effective lag (`Inf` for none).
#'
#' @return An object of class [sdvinecop-class].
#' @references Bladt, M. and McNeil, A. J. (2022). Time series copula models using d-vines
#'   and v-transforms. *Econometrics and Statistics*, 24, 27-48.
#'   \doi{10.1016/j.ecosta.2021.07.004}
#'
#'   Bladt, M. and McNeil, A. J. (2022). Time series with infinite-order partial
#'   copula dependence. *Dependence Modeling*, 10(1), 87-107.
#'   \doi{10.1515/demo-2022-0105}
#'
#'   Bladt, M., Dias, A., Han, J. and McNeil, A. J. (2026). Semiparametric and
#'   parametric distributional forecasting of univariate time series using
#'   non-Gaussian ARMA models based on D-vines. *Canadian Journal of Statistics*,
#'   54(1), e70033. \doi{10.1002/cjs.70033}
#'
#'   Dias, A., Han, J. and McNeil, A. J. (2027). GARCH copulas, v-transforms and
#'   D-vines for stochastic volatility. *Journal of Multivariate Analysis*, 217,
#'   105695. \doi{10.1016/j.jmva.2026.105695}
#'
#' @export
#'
#' @examples
#' # Gaussian ARMA copula with Gumbel at lag 1 and Clayton at lag 2
#' sdvinecop(armacop(ar = 0.7, ma = 0.3), pairs = list(pc("gumbel"), pc("clayton")))
#' # AR(2) with a t copula at lag 1: a finite-order d-vine
#' sdvinecop(armacop(ar = c(0.5, 0.2)), pairs = list(pc("t", aux = 6)))
sdvinecop <- function(kpacf = NULL, pairs = list(), tail = "gauss",
                     tautol = 1e-4, maxlag = Inf, vt = NULL, vtpar = "both") {
  if (!is.null(kpacf) && !is(kpacf, "gausscop")) {
    stop("'kpacf' must be a Gaussian copula process (an armacop or arfactcop object)")
  }
  if (is(pairs, "pc")) pairs <- list(pairs)
  if (!is.list(pairs)) stop("'pairs' must be a list")
  if (is(pairs, "bicop_dist")) pairs <- list(pairs)
  isbic <- vapply(pairs, is, logical(1), "bicop_dist")
  ispc <- vapply(pairs, is, logical(1), "pc")
  if (length(pairs) > 0L && !all(isbic | ispc)) {
    stop("'pairs' must contain pc() specifications or bicop_dist objects")
  }
  if (any(isbic) && any(ispc)) stop("cannot mix pc() specifications and bicop_dist objects in 'pairs'")
  free <- length(pairs) > 0L && all(isbic)
  if (free) {
    mode <- "free"
    if (!identical(tail, "gauss") && !(is(tail, "pc") && tail@family == "gauss")) {
      stop("free-mode models have a Gaussian tail; 'tail' cannot be changed")
    }
    tailspec <- NULL
  } else {
    mode <- "tied"
    if (is.null(kpacf)) stop("a tied model needs 'kpacf'")
    tailspec <- as_pc(tail)
  }
  K0 <- length(pairs)
  K <- NA_integer_
  limit <- "none"
  if (!is.null(kpacf) && gc_is_specified(kpacf)) {
    eff <- gc_effective_lag(kpacf, K0, tautol, maxlag)
    K <- eff$K
    limit <- eff$limit
  } else if (is.null(kpacf)) {
    K <- max(K0, 1L)
    limit <- "floor"
  }
  name <- if (mode == "tied") {
    sprintf("d-vine(%s, %d explicit, %s tail)", if (is.null(kpacf)) "" else kpacf@name,
            K0, tailspec@family)
  } else {
    sprintf("d-vine(free, %d explicit%s)", K0, if (is.null(kpacf)) "" else paste0(", ", kpacf@name, " tail"))
  }
  vtpar <- match.arg(vtpar, c("both", "common", "none"))
  has_ivt <- any(vapply(pairs, function(sp) is(sp, "pc") && sp@family == "ivt", logical(1))) ||
    (!is.null(tailspec) && tailspec@family == "ivt")
  if (has_ivt) {
    if (is.null(vt)) vt <- c(0.5, 0.5)
    if (!is.numeric(vt) || length(vt) != 2L) stop("'vt' must be two deltas in (0, 1)")
    vt <- as.numeric(vt)
    names(vt) <- c("delta1", "delta2")
  } else {
    if (!is.null(vt)) stop("'vt' is only used with inverse-v-transformed pair copulas (pc(\"ivt\", ...))")
    vt <- numeric(0)
  }
  new("sdvinecop", name = name, mode = mode, kpacf = kpacf, pairs = pairs, tail = tailspec,
      tautol = as.numeric(tautol), maxlag = as.numeric(maxlag), K = as.integer(K),
      Klimit = limit, vt = vt, vtpar = vtpar)
}

# Kendall's tau of the kpacf at lags 1..K for a specified Gaussian copula
# process.
dv_kpacf_tau <- function(G, K) gc_tau(G, K)

# Resolve a model to the list of K resolved pair copulas. Returns NULL if a
# pair copula specification is infeasible at the current parameters.
dv_resolve <- function(object) {
  K <- object@K
  if (is.na(K)) stop("the model has unspecified kpacf coefficients; fit it to data first")
  K0 <- length(object@pairs)
  pairs <- vector("list", K)
  if (object@mode == "tied") {
    tau <- dv_kpacf_tau(object@kpacf, K)
    vt <- if (length(object@vt)) object@vt else NULL
    for (j in seq_len(min(K0, K))) {
      pr <- dv_resolve_pair(object@pairs[[j]], tau[j], vt = vt)
      if (is.null(pr)) return(NULL)
      pairs[[j]] <- pr
    }
    if (K > K0) {
      tl <- dv_resolve_tail(object@tail, tau[(K0 + 1L):K], vt)
      if (is.null(tl)) return(NULL)
      pairs[(K0 + 1L):K] <- tl
    }
  } else {
    for (j in seq_len(min(K0, K))) pairs[[j]] <- dv_from_bicop(object@pairs[[j]])
    if (K > K0) {
      alpha <- gc_pacf(object@kpacf, K)
      for (j in (K0 + 1L):K) {
        pairs[[j]] <- list(family = "gauss", rotation = 0, parameters = alpha[j])
      }
    }
  }
  pairs
}

# Copula log-likelihood of complete data u in (0, 1) under a model; -Inf if the
# model is infeasible.
dv_model_loglik <- function(object, u) {
  pairs <- dv_resolve(object)
  if (is.null(pairs)) return(-Inf)
  dv_loglik(u, pairs)
}

#' @rdname tscop-methods
#' @export
setMethod("coef", "sdvinecop", function(object) {
  out <- if (is.null(object@kpacf)) numeric(0) else coef(object@kpacf)
  if (object@mode == "tied") {
    for (j in seq_along(object@pairs)) {
      sp <- object@pairs[[j]]
      if (sp@family %in% .pc_aux_families) {
        aux <- sp@aux
        names(aux) <- paste0(if (sp@family == "t") "nu" else "delta", j)
        out <- c(out, aux)
      }
    }
    if (length(object@vt) > 0L && object@vtpar != "none") {
      out <- c(out, if (object@vtpar == "common") c(delta = object@vt[[1L]]) else object@vt)
    }
  }
  out
})

# Description of the family used at each explicit lag, with the rotation
# resolved by the sign of tau when the ARMA coefficients are specified.
dv_explicit_labels <- function(object) {
  K0 <- length(object@pairs)
  if (K0 == 0L) return(character(0))
  if (object@mode == "free") {
    return(vapply(object@pairs, function(b) {
      rot <- if (b$rotation != 0) as.character(b$rotation) else ""
      paste0(b$family, rot)
    }, character(1)))
  }
  tau <- if (!is.null(object@kpacf) && gc_is_specified(object@kpacf)) dv_kpacf_tau(object@kpacf, K0) else rep(NA_real_, K0)
  vapply(seq_len(K0), function(j) {
    sp <- object@pairs[[j]]
    if (sp@family == "ivt") return(paste0("ivt(", sp@base, ")"))
    if (sp@family %in% .pc_radial) return(sp@family)
    if (is.na(tau[j])) return(paste0(sp@family, "(", sp@pos, "/", sp@neg, ")"))
    paste0(sp@family, if (tau[j] >= 0) sp@pos else sp@neg)
  }, character(1))
}

# Description of the tail family of a model: the base copula of an IVT tail and the rotations of the
# other non-radial families. Shared by the show methods of sdvinecop and sdvinecopfit.
dv_tail_label <- function(tail) {
  paste0(tail@family,
         if (tail@family == "ivt") paste0(" with base ", tail@base),
         if (!(tail@family %in% .pc_radial)) sprintf(" (positive: %g, negative: %g)", tail@pos, tail@neg))
}

#' @rdname tscop-methods
#' @export
setMethod("show", "sdvinecop", function(object) {
  cat("object class: sdvinecop (", object@mode, " mode)\n", sep = "")
  if (!is.null(object@kpacf)) cat("kpacf: ", gc_describe(object@kpacf), "\n", sep = "")
  lab <- dv_explicit_labels(object)
  if (length(lab) > 0L) cat(length(lab), " explicit pair copulas: ", paste(lab, collapse = ", "), "\n", sep = "")
  if (!is.null(object@tail)) cat("tail family: ", dv_tail_label(object@tail), "\n", sep = "")
  if (is.na(object@K)) {
    cat("effective lag: not yet determined (unspecified kpacf coefficients)\n")
  } else {
    basis <- switch(object@Klimit,
      ar = "pure AR order",
      tol = sprintf("tautol = %g", object@tautol),
      maxlag = sprintf("imposed maxlag = %g", object@maxlag),
      floor = "number of explicit pair copulas",
      cap = "internal limit reached (kpacf still above tautol)",
      none = "")
    cat("effective lag K = ", object@K, " (", basis, ")\n", sep = "")
  }
  if (length(object@vt) > 0L) {
    cat("v-transform deltas: ", paste(format(object@vt, digits = 4), collapse = ", "),
        switch(object@vtpar, both = "", common = " (common)", none = " (fixed)"), "\n", sep = "")
  }
  cf <- coef(object)
  if (length(cf) > 0L) {
    cat("parameters:\n")
    print(cf)
  }
  invisible(object)
})

# Complete copula data: strictly inside (0, 1), no missing values.
check_complete_copula_data <- function(U, min_n = 1L) {
  u <- check_copula_data(U, min_n = min_n)
  if (anyNA(u)) stop("d-vine models require complete data (no missing values)")
  u
}

# Resolved pair copulas of a model that must be fully specified and feasible.
dv_pairs_or_stop <- function(object) {
  pairs <- dv_resolve(object)
  if (is.null(pairs)) stop("the model is infeasible at its current parameter values")
  pairs
}

#' @describeIn tssim Exact stationary simulation from a known
#'   [sdvinecop-class] model, by the inverse Rosenblatt transform with the
#'   backward variables carried from step to step (cost O(nK)). `n` defaults
#'   to 1000.
#' @export
setMethod("tssim", "sdvinecop", function(object, n, ...) {
  if (missing(n)) n <- 1000
  if (!is.numeric(n) || length(n) != 1L || n < 1 || n != round(n)) {
    stop("'n' must be a positive whole number")
  }
  pairs <- dv_pairs_or_stop(object)
  as.numeric(dv_simulate(pairs, matrix(stats::runif(n), n, 1L)))
})

#' @describeIn ktau Kendall's tau of the pair copulas of a [sdvinecop-class]
#'   model at lags `1, ..., lagmax` (default 20); zero beyond the effective
#'   lag `K`. In tied mode the values are the ARMA kpacf.
#' @export
setMethod("ktau", "sdvinecop", function(object, lagmax, ...) {
  if (missing(lagmax)) lagmax <- 20
  if (!is.numeric(lagmax) || length(lagmax) != 1L || lagmax < 1 || lagmax != round(lagmax)) {
    stop("'lagmax' must be a positive whole number")
  }
  if (is.na(object@K)) stop("the model has unspecified kpacf coefficients; fit it to data first")
  K <- object@K
  tau <- numeric(lagmax)
  m <- min(K, lagmax)
  if (object@mode == "tied") {
    tau[seq_len(m)] <- dv_kpacf_tau(object@kpacf, K)[seq_len(m)]
  } else {
    pairs <- dv_pairs_or_stop(object)
    tau[seq_len(m)] <- vapply(pairs[seq_len(m)], dv_pair_ktau, numeric(1))
  }
  tau
})

dv_resid <- function(model, data, scale) {
  scale <- match.arg(scale, c("uniform", "normal"))
  u <- check_complete_copula_data(data)
  r <- dv_rosenblatt(dv_pairs_or_stop(model), u)
  out <- if (scale == "uniform") r else stats::qnorm(r)
  restore_attrs(out, data)
}

#' @rdname resid-method
#' @export
setMethod("resid", "sdvinecop", function(object, data, scale = "uniform", ...) {
  dv_resid(object, data, scale)
})

#' @rdname resid-method
#' @export
setMethod("residuals", "sdvinecop", function(object, data, scale = "uniform", ...) {
  dv_resid(object, data, scale)
})

dv_predict <- function(model, data, h, type, x) {
  type <- match.arg(type, c("df", "qf", "dens"))
  if (!identical(as.numeric(h), 1)) stop("only h = 1 is implemented for d-vine models")
  if (is.null(x) || !is.numeric(x) || anyNA(x) || any(x <= 0 | x >= 1)) {
    stop("'x' must be numeric values strictly inside (0, 1)")
  }
  u <- check_complete_copula_data(data)
  dv_predictive(dv_pairs_or_stop(model), u, as.numeric(x), type)
}

#' @rdname predict-method
#' @export
setMethod("predict", "sdvinecop", function(object, data, h = 1, type = "df", x = NULL, ...) {
  dv_predict(object, data, h, type, x)
})
