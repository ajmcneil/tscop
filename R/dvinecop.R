#' Stationary d-vine copula processes
#'
#' Stationary d-vine copula process of finite order \eqn{K}, built from a
#' sequence of pair copulas \eqn{C_1, \ldots, C_K} (independence beyond lag
#' \eqn{K}).
#'
#' There are two modes.
#'
#' **Tied mode** (the only mode that can be fitted). A known Gaussian
#' ARMA/SARMA [armacop-class] model supplies the Kendall partial
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
#' @slot arma an [armacop-class] object, or `NULL`.
#' @slot pairs list of explicit pair copulas: [pc-class] objects (tied mode)
#'   or `bicop_dist` objects (free mode).
#' @slot tail the tail family: a [pc-class] object (tied mode), or `NULL`.
#' @slot tautol,maxlag the tolerance and cap defining the effective lag.
#' @slot K the effective lag (`NA` for a model structure whose ARMA
#'   coefficients are not yet specified).
#' @slot Klimit what determined `K`: `"ar"` (pure AR order), `"tol"`,
#'   `"maxlag"`, `"floor"` (number of explicit pair copulas), `"cap"`
#'   (internal limit on the lags examined) or `"none"` (free mode without an
#'   ARMA tail).
#'
#' @exportClass dvinecop
setClass("dvinecop", slots = c(
  name = "character", mode = "character", arma = "armacopOrNULL",
  pairs = "list", tail = "pcOrNULL", tautol = "numeric", maxlag = "numeric",
  K = "integer", Klimit = "character"
))

setValidity("dvinecop", function(object) {
  if (!(object@mode %in% c("tied", "free"))) return("'mode' must be \"tied\" or \"free\"")
  if (length(object@tautol) != 1L || !(object@tautol > 0)) return("'tautol' must be a positive number")
  if (length(object@maxlag) != 1L || is.na(object@maxlag) || object@maxlag < 1) {
    return("'maxlag' must be at least 1 (use Inf for no cap)")
  }
  K0 <- length(object@pairs)
  if (object@mode == "tied") {
    if (is.null(object@arma)) return("a tied model needs an 'arma' component")
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
    if (K0 == 0L && is.null(object@arma)) return("specify 'pairs' or 'arma'")
  }
  if (!is.na(object@K) && object@K < max(K0, 1L)) return("'K' must be at least the number of explicit pairs")
  TRUE
})

#' Constructor for stationary d-vine copula processes
#'
#' Creates a [dvinecop-class] object in tied mode (`pairs` are
#' [pc()] specifications) or free mode (`pairs` are
#' [rvinecopulib::bicop_dist()] objects).
#'
#' @param arma an [armacop-class] object. In tied mode its Kendall partial
#'   autocorrelations calibrate every pair copula; it may be a model structure
#'   (`armacop(order = ...)`) if the model is only to be passed to a fitting
#'   function. In free mode it is optional and supplies a Gaussian tail.
#' @param pairs list of explicit pair copulas for lags `1, ..., length(pairs)`:
#'   all [pc()] specifications (tied mode) or all `bicop_dist` objects (free
#'   mode). A `pc("gauss")` entry keeps a lag Gaussian.
#' @param tail the family used for all lags beyond `pairs` in tied mode:
#'   `"gauss"` (default), `"frank"`, `"clayton"`, `"gumbel"` or `"joe"`, or a
#'   [pc()] specification of one of these. Families whose Kendall's tau cannot
#'   tend to zero with their auxiliary parameters fixed (`"t"`, `"bb1"`) are
#'   rejected. Free mode has a Gaussian tail.
#' @param tautol tolerance on the absolute Kendall's tau defining the
#'   effective lag.
#' @param maxlag cap on the effective lag (`Inf` for none).
#'
#' @return An object of class [dvinecop-class].
#' @export
#'
#' @examples
#' # Gaussian ARMA copula with Gumbel at lag 1 and Clayton at lag 2
#' dvinecop(armacop(ar = 0.7, ma = 0.3), pairs = list(pc("gumbel"), pc("clayton")))
#' # AR(2) with a t copula at lag 1: a finite-order d-vine
#' dvinecop(armacop(ar = c(0.5, 0.2)), pairs = list(pc("t", aux = 6)))
dvinecop <- function(arma = NULL, pairs = list(), tail = "gauss",
                     tautol = 1e-4, maxlag = Inf) {
  if (!is.null(arma) && !is(arma, "armacop")) stop("'arma' must be an armacop object")
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
    if (is.null(arma)) stop("a tied model needs 'arma'")
    tailspec <- as_pc(tail)
  }
  K0 <- length(pairs)
  K <- NA_integer_
  limit <- "none"
  if (!is.null(arma) && arma_is_specified(arma)) {
    e <- arma_coef_expanded(arma)
    eff <- dv_effective_lag(e$ar, e$ma, K0, tautol, maxlag, arma@period)
    K <- eff$K
    limit <- eff$limit
  } else if (is.null(arma)) {
    K <- max(K0, 1L)
    limit <- "floor"
  }
  name <- if (mode == "tied") {
    sprintf("d-vine(%s, %d explicit, %s tail)", if (is.null(arma)) "" else arma@name,
            K0, tailspec@family)
  } else {
    sprintf("d-vine(free, %d explicit%s)", K0, if (is.null(arma)) "" else paste0(", ", arma@name, " tail"))
  }
  new("dvinecop", name = name, mode = mode, arma = arma, pairs = pairs, tail = tailspec,
      tautol = as.numeric(tautol), maxlag = as.numeric(maxlag), K = as.integer(K),
      Klimit = limit)
}

# Kendall's tau of the ARMA kpacf at lags 1..K for a specified armacop.
dv_arma_tau <- function(arma, K) {
  e <- arma_coef_expanded(arma)
  (2 / pi) * asin(arma_pacf(e$ar, e$ma, K))
}

# Resolve a model to the list of K resolved pair copulas. Returns NULL if a
# pair copula specification is infeasible at the current parameters.
dv_resolve <- function(object) {
  K <- object@K
  if (is.na(K)) stop("the model has unspecified ARMA coefficients; fit it to data first")
  K0 <- length(object@pairs)
  pairs <- vector("list", K)
  if (object@mode == "tied") {
    tau <- dv_arma_tau(object@arma, K)
    for (j in seq_len(K)) {
      spec <- if (j <= K0) object@pairs[[j]] else object@tail
      pr <- dv_resolve_pair(spec, tau[j])
      if (is.null(pr)) return(NULL)
      pairs[[j]] <- pr
    }
  } else {
    for (j in seq_len(min(K0, K))) pairs[[j]] <- dv_from_bicop(object@pairs[[j]])
    if (K > K0) {
      e <- arma_coef_expanded(object@arma)
      alpha <- arma_pacf(e$ar, e$ma, K)
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

#' @describeIn dvinecop-class Coefficients: the ARMA coefficients followed by
#'   the auxiliary parameters (`nu<lag>` for t copulas, `delta<lag>` for BB1
#'   copulas) of the explicit specifications (tied mode).
#' @param object a [dvinecop-class] object.
#' @export
setMethod("coef", "dvinecop", function(object) {
  out <- if (is.null(object@arma)) numeric(0) else coef(object@arma)
  if (object@mode == "tied") {
    for (j in seq_along(object@pairs)) {
      sp <- object@pairs[[j]]
      if (sp@family %in% .pc_aux_families) {
        aux <- sp@aux
        names(aux) <- paste0(if (sp@family == "t") "nu" else "delta", j)
        out <- c(out, aux)
      }
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
  tau <- if (!is.null(object@arma) && arma_is_specified(object@arma)) dv_arma_tau(object@arma, K0) else rep(NA_real_, K0)
  vapply(seq_len(K0), function(j) {
    sp <- object@pairs[[j]]
    if (sp@family %in% .pc_radial) return(sp@family)
    if (is.na(tau[j])) return(paste0(sp@family, "(", sp@pos, "/", sp@neg, ")"))
    paste0(sp@family, if (tau[j] >= 0) sp@pos else sp@neg)
  }, character(1))
}

#' @describeIn dvinecop-class Show method.
#' @export
setMethod("show", "dvinecop", function(object) {
  cat("object class: dvinecop (", object@mode, " mode)\n", sep = "")
  if (!is.null(object@arma)) {
    cat("ARMA kpacf: ", object@arma@name,
        if (!is.na(object@arma@period)) paste0(" with period ", object@arma@period), "\n", sep = "")
  }
  lab <- dv_explicit_labels(object)
  if (length(lab) > 0L) cat(length(lab), " explicit pair copulas: ", paste(lab, collapse = ", "), "\n", sep = "")
  if (!is.null(object@tail)) cat("tail family: ", object@tail@family,
                                 if (!(object@tail@family %in% .pc_radial))
                                   sprintf(" (positive: %g, negative: %g)", object@tail@pos, object@tail@neg),
                                 "\n", sep = "")
  if (is.na(object@K)) {
    cat("effective lag: not yet determined (unspecified ARMA coefficients)\n")
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
  cf <- coef(object)
  if (length(cf) > 0L) {
    cat("parameters:\n")
    print(cf)
  }
  invisible(object)
})
