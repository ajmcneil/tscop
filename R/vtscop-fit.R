# Maximum likelihood fitting of vtscop models.

#' Fitted process-level v-transform copula processes
#'
#' Class of objects returned by [tsfit()] for [vtscop-class] models.
#'
#' @slot model the fitted [vtscop-class] model.
#' @slot data the copula data to which the model was fitted.
#' @slot pseudo logical: were the data pseudo-observations?
#' @slot estimate named numeric vector of estimates (base coefficients, then
#'   the free v-transform parameters).
#' @slot se standard errors matching `estimate`; `NA` when not computed.
#' @slot vcov the matching covariance matrix; `NA` when not computed.
#' @slot se_method `"none"`, `"hessian"` or `"bootstrap"`.
#' @slot boot matrix of bootstrap estimates (no rows unless
#'   `se = "bootstrap"`).
#' @slot loglik maximised copula log-likelihood.
#' @slot nobs number of non-missing observations.
#' @slot npar number of estimated parameters.
#' @slot convergence convergence code of the last optimisation (0 is success).
#' @slot message optimiser message.
#' @slot start the [vtscop-class] model with the starting values.
#' @slot Kstart,nrefit,dloglik effective-lag bookkeeping for a [sdvinecop-class]
#'   base (as in [sdvinecopfit-class]); `NA` for an `arfactcop` base.
#' @slot template the original input if it was a `ts`, `zoo` or `xts` object,
#'   otherwise `NULL`.
#'
#' @exportClass vtscopfit
setClass("vtscopfit", slots = c(
  model = "vtscop", data = "numeric", pseudo = "logical",
  estimate = "numeric", se = "numeric", vcov = "matrix", se_method = "character",
  boot = "matrix", loglik = "numeric", nobs = "integer", npar = "integer",
  convergence = "integer", message = "character", start = "vtscop",
  Kstart = "integer", nrefit = "integer", dloglik = "numeric", template = "ANY"
))

# ---- layout and objective ---------------------------------------------------

vt_layout <- function(object) {
  base <- object@base
  p <- vt_params(object@vt)
  free <- setdiff(names(p), object@vtfix)
  if (is(base, "arfactcop")) {
    ops <- gc_free_ops(base)
    list(kind = "afc", ops = ops, lay = NULL, nb = ops$nfree, m = base@m, noise = base@noise,
         vtname = object@vt@name, vt_pars = p, free = free)
  } else {
    lay <- dv_layout(base)
    list(kind = "dvine", ops = NULL, lay = lay, nb = lay$nk + length(lay$aux_lags),
         vtname = object@vt@name, vt_pars = p, free = free)
  }
}

# Full v-transform parameters from the free-scale values: fulcrum on the
# logit scale, kappa and xi on the log scale.
vt_pars_from_free <- function(zv, L) {
  p <- L$vt_pars
  for (i in seq_along(L$free)) {
    nm <- L$free[i]
    p[[nm]] <- if (nm == "delta") stats::plogis(zv[i]) else exp(zv[i])
  }
  p
}

vt_pars_to_free <- function(p, free) {
  vapply(free, function(nm) if (nm == "delta") stats::qlogis(min(max(p[[nm]], 1e-3), 1 - 1e-3)) else log(p[[nm]]),
         numeric(1), USE.NAMES = FALSE)
}

# Negative copula log-likelihood on the free scale (base free vector, then the
# v-transform free parameters) at fixed effective lag K (sdvinecop bases).
vt_negll <- function(z, L, K, u) {
  nb <- L$nb
  zv <- z[nb + seq_along(L$free)]
  if (any(!is.finite(z)) || (length(zv) > 0L && any(abs(zv) > 8))) return(1e10)
  p <- vt_pars_from_free(zv, L)
  v <- tryCatch(vt_transform_data(vt_build(L$vtname, p), u), error = function(e) NULL)
  if (is.null(v) || anyNA(v[!is.na(u)])) return(1e10)
  zb <- z[seq_len(nb)]
  if (L$kind == "afc") {
    x <- stats::qnorm(v)
    val <- afc_negll(zb, L$m, L$noise, x)
    if (val >= 1e10) return(1e10)
    val + sum(stats::dnorm(x, log = TRUE), na.rm = TRUE)
  } else {
    dv_negll(zb, L$lay, K, v)
  }
}

# Gradient of vt_negll for a d-vine base from the C++ kernel. The base parameters
# enter through the level parameters (as in dv_negll_grad); the v-transform
# parameters enter through the transformed series v = V(u), whose derivatives
# with respect to the free v-transform parameters (central differences of the
# udp transform, chained through the logit and log maps) are passed to the
# kernel as derivatives of the data. NULL if the kernel does not apply (a
# Gaussian base, unsupported pair copulas, more than 12 free parameters).
vt_negll_grad <- function(z, L, K, u, du = NULL) {
  if (L$kind == "afc") return(vt_negll_grad_afc(z, L, u, du))
  if (L$kind != "dvine" || anyNA(K)) return(NULL)
  nb <- L$nb
  nv <- length(L$free)
  zv <- z[nb + seq_len(nv)]
  nx <- if (is.null(du)) 0L else ncol(du)    # extra leading directions: derivatives of the data u
  if (any(!is.finite(z)) || (nv > 0L && any(abs(zv) > 8)) || nx + nb + nv > 12L) return(NULL)
  g <- dv_grad_pieces(z[seq_len(nb)], L$lay, K)
  if (is.null(g)) return(NULL)
  p <- vt_pars_from_free(zv, L)
  vt <- vt_build(L$vtname, p)
  v <- vt_transform_data(vt, u)
  if (anyNA(v)) return(NULL)
  P <- nx + nb + nv
  dvm <- matrix(0, length(u), P)
  if (nx > 0L) {
    # v depends on u through the transform: dv/du by central differences, times du/d(extra directions)
    h <- 1e-6
    up <- pmin(u + h, 1 - 1e-12)
    um <- pmax(u - h, 1e-12)
    dvdu <- (vt_transform_data(vt, up) - vt_transform_data(vt, um)) / (up - um)
    dvm[, seq_len(nx)] <- dvdu * du
  }
  for (k in seq_len(nv)) {
    nm <- L$free[k]
    th <- p[[nm]]
    h <- 1e-6
    pp <- p; pp[[nm]] <- th + h
    pm <- p; pm[[nm]] <- th - h
    dv <- (vt_transform_data(vt_build(L$vtname, pp), u) - vt_transform_data(vt_build(L$vtname, pm), u)) / (2 * h)
    dvm[, nx + nb + k] <- dv * (if (nm == "delta") th * (1 - th) else th)
  }
  e <- g$e
  pad <- function(m) cbind(matrix(0, nrow(m), nx), m, matrix(0, nrow(m), nv))
  r <- dv_loglik_grad_cpp(v, dvm, e$fam, e$rot, e$par, pad(g$dpar), e$ivt, e$bfam, e$brot, e$bpar, pad(g$dbpar),
                          e$d1, e$d2, c(numeric(nx), g$dd1, numeric(nv)), c(numeric(nx), g$dd2, numeric(nv)))
  if (!is.finite(r$value) || any(!is.finite(r$gradient))) return(NULL)
  -r$gradient
}

# The same for an arfactcop base: the base is a Gaussian copula, so the objective is
# the Kalman log-likelihood of x = qnorm(v) plus the sum of log phi(x); the
# derivatives of x along the free v-transform parameters (and the extra data
# directions du) go through the filter as derivatives of the data.
vt_negll_grad_afc <- function(z, L, u, du = NULL) {
  if (anyNA(u)) return(NULL)
  nb <- L$nb
  nv <- length(L$free)
  nx <- if (is.null(du)) 0L else ncol(du)
  zv <- z[nb + seq_len(nv)]
  if (any(!is.finite(z)) || (nv > 0L && any(abs(zv) > 8))) return(NULL)
  p <- vt_pars_from_free(zv, L)
  vt <- vt_build(L$vtname, p)
  v <- vt_transform_data(vt, u)
  if (anyNA(v)) return(NULL)
  x <- stats::qnorm(v)
  P <- nx + nb + nv
  dvm <- matrix(0, length(u), P)
  if (nx > 0L) {
    h <- 1e-6
    up <- pmin(u + h, 1 - 1e-12)
    um <- pmax(u - h, 1e-12)
    dvm[, seq_len(nx)] <- (vt_transform_data(vt, up) - vt_transform_data(vt, um)) / (up - um) * du
  }
  for (k in seq_len(nv)) {
    nm <- L$free[k]
    th <- p[[nm]]
    h <- 1e-6
    pp <- p; pp[[nm]] <- th + h
    pm <- p; pm[[nm]] <- th - h
    dv <- (vt_transform_data(vt_build(L$vtname, pp), u) - vt_transform_data(vt_build(L$vtname, pm), u)) / (2 * h)
    dvm[, nx + nb + k] <- dv * (if (nm == "delta") th * (1 - th) else th)
  }
  dx <- dvm / stats::dnorm(x)
  g <- afc_negll_grad(z[seq_len(nb)], L$m, L$noise, x, dx, before = nx, after = nv)
  if (is.null(g)) return(NULL)
  g - colSums(x * dx)
}

# Gradient for optim(): exact where available, otherwise central differences of vt_negll.
vt_negll_gr <- function(z, L, K, u) {
  g <- vt_negll_grad(z, L, K, u)
  if (!is.null(g)) return(g)
  vapply(seq_along(z), function(k) {
    h <- 1e-6
    zp <- z; zp[k] <- zp[k] + h
    zm <- z; zm[k] <- zm[k] - h
    (vt_negll(zp, L, K, u) - vt_negll(zm, L, K, u)) / (2 * h)
  }, numeric(1))
}

# Natural coefficient vector (base coefficients, then the free v-transform
# parameters) from the free vector.
vt_free_to_coef <- function(z, L) {
  nb <- L$nb
  zb <- z[seq_len(nb)]
  base <- if (L$kind == "afc") L$ops$coef(zb) else dv_free_to_coef(zb, L$lay)
  p <- vt_pars_from_free(z[nb + seq_along(L$free)], L)
  c(base, unlist(p[L$free]))
}

# ---- starting values -------------------------------------------------------

vt_start <- function(object, init, u, control) {
  L <- vt_layout(object)
  nb <- L$nb
  p0 <- L$vt_pars
  base <- object@base
  ok <- !is.na(u)
  if (!is.null(init)) {
    if (is(init, "vtscop")) {
      if (!identical(class(init@base), class(base)) || !identical(init@vt@name, object@vt@name)) {
        stop("'init' must be a vtscop with the same base class and v-transform as 'object'")
      }
      p0 <- vt_params(init@vt)
      init_base <- init@base
      init_coef <- NULL
    } else if (is.numeric(init) && length(init) == length(coef(base)) + length(L$free)) {
      nbc <- length(coef(base))
      p0[L$free] <- init[nbc + seq_along(L$free)]
      init_coef <- as.numeric(init[seq_len(nbc)])
      init_base <- NULL
    } else {
      stop("'init' must be a vtscop model or a numeric vector of coefficients")
    }
  } else {
    init_base <- NULL
    init_coef <- NULL
  }
  v0 <- rep(NA_real_, length(u))
  v0[ok] <- vt_transform_data(vt_build(L$vtname, p0), u[ok])
  if (L$kind == "afc") {
    kp <- if (!is.null(init_base)) init_base@base else if (!is.null(init_coef)) gc_from_coef(base, init_coef) else base
    if (!gc_is_specified(kp)) kp <- gc_prefit(kp, v0, control)
    zb <- L$ops$to_free(kp)
    st <- NULL
  } else {
    st <- dv_start(base, if (!is.null(init_base)) init_base@base else init_coef, v0, control)
    zb <- st$z0
  }
  zv <- if (length(L$free) > 0L) vt_pars_to_free(p0, L$free) else numeric(0)
  list(z0 = c(zb, zv), L = L, vt0 = p0, st = st)
}

# ---- fitting procedure --------------------------------------------------------

vt_fit_all <- function(u, object, init, control) {
  stt <- vt_start(object, init, u, control)
  L <- stt$L
  z0 <- stt$z0
  ctrl <- modifyList(list(maxit = 500, reltol = 1e-10), control)
  fit_fixed <- function(K, z) {
    npar <- length(z)
    if (npar == 0L) return(list(z = z, value = vt_negll(z, L, K, u), convergence = 0L, message = ""))
    c1 <- modifyList(ctrl, list(ndeps = rep(1e-6, npar)))
    opt <- stats::optim(z, vt_negll, gr = vt_negll_gr, L = L, K = K, u = u, method = "BFGS", control = c1)
    if (opt$convergence != 0L || opt$value >= 1e10) {
      alt <- stats::optim(opt$par, vt_negll, L = L, K = K, u = u, method = "Nelder-Mead",
                          control = list(maxit = 5000, reltol = 1e-10, warn.1d.NelderMead = FALSE))
      if (alt$value < opt$value) opt <- alt
    }
    list(z = opt$par, value = opt$value, convergence = opt$convergence,
         message = if (is.null(opt$message)) "" else opt$message)
  }
  if (L$kind == "afc") {
    if (vt_negll(z0, L, NA, u) >= 1e10) stop("the starting values are infeasible for these data")
    fit <- fit_fixed(NA, z0)
    return(list(z = fit$z, K = NA_integer_, Kstart = NA_integer_, nrefit = NA_integer_, dloglik = NA_real_,
                limit = "none", fit = fit, L = L, stt = stt, loglik = -vt_negll(fit$z, L, NA, u)))
  }
  base <- object@base
  lay <- L$lay
  tautol <- base@tautol
  maxlag <- base@maxlag
  zb_idx <- seq_len(L$nb)
  K_start <- dv_K_at(z0[zb_idx], lay, tautol, maxlag)$K
  if (vt_negll(z0, L, K_start, u) >= 1e10) stop("the starting values are infeasible for these data")
  K_of <- function(z) dv_K_at(z[zb_idx], lay, tautol, maxlag)$K
  ll_of <- function(z, K) -vt_negll(z, L, K, u)
  res <- dv_k_rule(fit_fixed, K_of, ll_of, K_start, z0)
  limit <- dv_K_at(res$z[zb_idx], lay, tautol, maxlag)$limit
  if (res$K != K_of(res$z)) limit <- "tol"
  c(res, list(L = L, stt = stt, limit = limit, loglik = ll_of(res$z, res$K)))
}

# The fitted model from a free vector.
vt_build_model <- function(object, z, L, K, limit) {
  nb <- L$nb
  zb <- z[seq_len(nb)]
  base <- if (L$kind == "afc") L$ops$build(zb) else dv_build_model(object@base, zb, L$lay, K, limit)
  p <- vt_pars_from_free(z[nb + seq_along(L$free)], L)
  vtscop(base, vt_build(L$vtname, p), vtfix = object@vtfix)
}

#' @describeIn tsfit Maximum likelihood fit of a [vtscop-class] model: the base
#'   parameters and the free v-transform parameters (fulcrum, kappa, xi) are
#'   estimated jointly on an unconstrained scale (BFGS with a Nelder-Mead
#'   fallback). The base is first fitted as a Gaussian copula (for a d-vine base
#'   with the effective-lag rule of [sdvinecop-class]) to the v-transformed
#'   data at the starting v-transform. `init` is a [vtscop-class] model of the
#'   same structure or a numeric vector in `coef` order; `bootstart` is
#'   accepted for compatibility: bootstrap refits start from the fitted values
#'   (with the effective lag fixed for a d-vine base).
#'
#' @return For a [vtscop-class] model, an object of class [vtscopfit-class].
#' @export
setMethod("tsfit", signature(U = "ANY", object = "vtscop"),
  function(U, object, pseudo = TRUE, se = FALSE, B = 200, init = NULL,
           control = list(), bootstart = "fitted", cores = 1L, ...) {
    if (isTRUE(se)) se <- "bootstrap"
    if (!identical(se, FALSE) && !(is.character(se) && length(se) == 1L &&
                                   se %in% c("hessian", "bootstrap"))) {
      stop("'se' must be FALSE, \"hessian\" or \"bootstrap\"")
    }
    if (identical(se, "bootstrap") &&
        (!is.numeric(B) || length(B) != 1L || B < 10 || B != round(B))) {
      stop("'B' must be a whole number of at least 10")
    }
    u <- if (is(object@base, "sdvinecop")) check_complete_copula_data(U, min_n = 20L) else check_copula_data(U)
    n <- length(u)
    fit <- vt_fit_all(u, object, init, control)
    L <- fit$L
    z <- fit$z
    K <- fit$K
    model <- vt_build_model(object, z, L, K, fit$limit)
    est <- coef(model)
    cn <- names(est)
    npar <- length(est)
    V <- matrix(NA_real_, npar, npar, dimnames = list(cn, cn))
    boot <- matrix(NA_real_, 0L, npar, dimnames = list(NULL, cn))
    se_method <- "none"
    if (identical(se, "hessian") && npar > 0L) {
      H <- num_hessian(function(v) vt_negll(v, L, K, u), z)
      Vz <- tryCatch(solve(H), error = function(e) NULL)
      if (is.null(Vz) || any(!is.finite(Vz)) || any(diag(Vz) <= 0) ||
          any(eigen(H, symmetric = TRUE, only.values = TRUE)$values <= 0)) {
        warning("the information matrix is not positive definite; standard errors not available")
      } else {
        J <- num_jacobian(function(v) vt_free_to_coef(v, L), z)
        Vn <- J %*% Vz %*% t(J)
        V[] <- (Vn + t(Vn)) / 2
        se_method <- "hessian"
      }
    } else if (identical(se, "bootstrap") && npar > 0L) {
      Ub <- vt_simulate(model, n, B)
      boot <- matrix(NA_real_, B, npar, dimnames = list(NULL, cn))
      ctrl <- modifyList(list(maxit = 500, reltol = 1e-10, ndeps = rep(1e-6, length(z))), control)
      miss <- is.na(u)
      Ub[miss, ] <- NA
      if (pseudo) Ub <- apply(Ub, 2L, function(v) rank(v, na.last = "keep", ties.method = "random") / (sum(!miss) + 1))
      Hinv <- boot_prepare(function(zz) vt_negll_gr(zz, L, K, u), z)
      refit <- function(b) {
        ub <- Ub[, b]
        tryCatch(suppressWarnings({
          r <- if (is.null(Hinv)) NULL else boot_newton(z, function(zz) vt_negll(zz, L, K, ub),
                                                        function(zz) vt_negll_gr(zz, L, K, ub), Hinv)
          if (is.null(r)) {
            opt <- stats::optim(z, vt_negll, gr = vt_negll_gr, L = L, K = K, u = ub, method = "BFGS", control = ctrl)
            if (opt$convergence != 0L || opt$value >= 1e10) stop("no convergence")
            r <- list(par = opt$par)
          }
          vt_free_to_coef(r$par, L)
        }), error = function(e) NULL)
      }
      res <- boot_run(B, refit, cores)
      for (b in seq_len(B)) if (is.numeric(res[[b]]) && length(res[[b]]) == npar) boot[b, ] <- res[[b]]
      good <- stats::complete.cases(boot)
      if (sum(good) >= 2L) {
        V[] <- stats::cov(boot[good, , drop = FALSE])
        se_method <- "bootstrap"
        if (sum(good) < B) warning(sprintf("%d of %d bootstrap refits failed", B - sum(good), B))
      } else {
        warning("too few successful bootstrap refits; standard errors not available")
      }
    }
    se_vec <- if (se_method == "none") rep(NA_real_, npar) else sqrt(diag(V))
    names(se_vec) <- cn
    start_model <- tryCatch(vt_build_model(object, fit$stt$z0, L, if (is.na(fit$Kstart)) NA_integer_ else fit$Kstart, "tol"),
                            error = function(e) object)
    new("vtscopfit",
      model = model, data = u, pseudo = pseudo, estimate = est, se = se_vec, vcov = V,
      se_method = se_method, boot = boot, loglik = fit$loglik, nobs = as.integer(sum(!is.na(u))),
      npar = as.integer(npar), convergence = as.integer(fit$fit$convergence),
      message = fit$fit$message, start = start_model,
      Kstart = as.integer(fit$Kstart), nrefit = as.integer(fit$nrefit), dloglik = as.numeric(fit$dloglik),
      template = if (inherits(U, c("ts", "zoo", "xts"))) U else NULL
    )
  }
)

# ---- methods for fitted objects -----------------------------------------------

#' @rdname tscop-methods
#' @export
setMethod("coef", "vtscopfit", function(object) object@estimate)

#' @rdname tscop-methods
#' @export
setMethod("logLik", "vtscopfit", function(object) {
  structure(object@loglik, df = object@npar, nobs = object@nobs, class = "logLik")
})

#' @rdname tscop-methods
#' @export
setMethod("vcov", "vtscopfit", function(object, ...) object@vcov)

#' @rdname tscop-methods
#' @export
setMethod("show", "vtscopfit", function(object) {
  m <- object@model
  cat("object class: vtscopfit\n")
  p <- vt_params(m@vt)
  cat("v-transform: ", m@vt@name, if (length(p) > 0L) paste0(" (", paste(names(p), "=", format(p, digits = 4), collapse = ", "), ")"), "\n", sep = "")
  b <- utils::capture.output(show(m@base))
  cat("base process:\n"); cat(paste0("  ", b), sep = "\n")
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
  invisible(object)
})

#' @rdname tssim
#' @export
setMethod("tssim", "vtscopfit", function(object, n, ...) {
  if (missing(n)) n <- length(object@data)
  tssim(object@model, n)
})

#' @rdname ktau
#' @export
setMethod("ktau", "vtscopfit", function(object, lagmax, ...) {
  if (missing(lagmax)) lagmax <- 20
  ktau(object@model, lagmax)
})
