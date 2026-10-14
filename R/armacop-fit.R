# Maximum likelihood fitting of armacop models.

#' Fitted ARMA copula processes
#'
#' Class of objects returned by [tsfit()] for [armacop-class] models.
#'
#' @slot model the fitted [armacop-class] model.
#' @slot data the copula data to which the model was fitted.
#' @slot pseudo logical: were the data pseudo-observations?
#' @slot estimate named numeric vector of coefficient estimates.
#' @slot se standard errors matching `estimate`; `NA` when not computed.
#' @slot vcov the matching covariance matrix; `NA` when not computed.
#' @slot se_method `"none"`, `"hessian"` or `"bootstrap"`.
#' @slot boot matrix of bootstrap estimates (no rows unless
#'   `se = "bootstrap"`).
#' @slot loglik maximised copula log-likelihood.
#' @slot nobs number of non-missing observations.
#' @slot npar number of estimated coefficients.
#' @slot convergence convergence code of the optimiser (0 is success).
#' @slot message optimiser message.
#' @slot start the [armacop-class] model used as starting values.
#' @slot template the original input if it was a `ts`, `zoo` or `xts` object
#'   (used to restore the time index on residuals), otherwise `NULL`.
#'
#' @exportClass armacopfit
setClass("armacopfit", slots = c(
  model = "armacop", data = "numeric", pseudo = "logical",
  estimate = "numeric", se = "numeric", vcov = "matrix", se_method = "character",
  boot = "matrix", loglik = "numeric", nobs = "integer", npar = "integer",
  convergence = "integer", message = "character", start = "armacop",
  template = "ANY"
))

# Validate copula data: numeric, in (0, 1), NA allowed.
check_copula_data <- function(U, min_n = 20L) {
  if (!is.numeric(U)) stop("'U' must be numeric")
  u <- as.numeric(U)
  ok <- !is.na(u)
  if (any(!is.finite(u[ok]))) stop("'U' must be finite")
  if (any(u[ok] <= 0 | u[ok] >= 1)) {
    stop("'U' must lie strictly inside (0, 1); see udp::boundaryadjust() for data on the boundary")
  }
  if (sum(ok) < min_n) {
    stop(sprintf("at least %d non-missing observation%s required", min_n,
                 if (min_n == 1L) " is" else "s are"))
  }
  u
}

# Split a free (unconstrained) vector into the four coefficient blocks.
split_blocks <- function(v, ord) {
  idx <- rep(c("ar", "ma", "sar", "sma"), ord)
  split(v, factor(idx, levels = c("ar", "ma", "sar", "sma")))
}

# Unconstrained vector -> expanded ARMA coefficients.
free_to_expanded <- function(z, ord, period) {
  b <- split_blocks(z, ord)
  arma_expand(
    ar = free_to_ar(b$ar), ma = free_to_ma(b$ma),
    sar = free_to_ar(b$sar), sma = free_to_ma(b$sma), period = period
  )
}

# Natural coefficient list -> unconstrained vector.
pars_to_free <- function(pars) {
  c(ar_to_free(pars$ar), ma_to_free(pars$ma),
    ar_to_free(pars$sar), ma_to_free(pars$sma))
}

# Negative Gaussian log-likelihood on the unconstrained scale. The Jacobian
# term converting it to a copula likelihood does not depend on the parameters.
armacop_negll <- function(z, ord, period, x) {
  if (any(!is.finite(z)) || max(abs(z)) > 6 + 1e-8) return(1e10)
  e <- free_to_expanded(z, ord, period)
  val <- tryCatch(-arma_loglik(x, e$ar, e$ma), error = function(err) NA_real_)
  if (!is.finite(val)) 1e10 else val
}

# Gaussian log-likelihood of the normal scores: the compiled filter
# (src/arma_kalman.cpp) unless options(tscop.cpp = FALSE), else the R filter.
arma_loglik <- function(x, ar, ma, tol = 1e-12) {
  if (isTRUE(getOption("tscop.cpp", TRUE))) arma_loglik_cpp(as.numeric(x), ar, ma, tol)
  else arma_filter(x, ar, ma, tol)$loglik
}

# Gradient of armacop_negll from the compiled filter: the derivatives of the expanded
# AR and MA coefficients with respect to the free vector z by central differences (the
# map is smooth), the kernel gradient along them. The directions are `before` extra
# data directions, then z; dx is the n x (before + length(z)) matrix of the derivatives
# of x along them (NULL for constant data). NULL if unavailable.
armacop_negll_grad <- function(z, ord, period, x, dx = NULL, before = 0L) {
  if (!isTRUE(getOption("tscop.cpp", TRUE)) || !isTRUE(getOption("tscop.cppgrad", TRUE))) return(NULL)
  nz <- length(z)
  P <- before + nz
  if (nz < 1L || P > 12L || any(!is.finite(z)) || max(abs(z)) > 6 + 1e-8) return(NULL)
  e <- free_to_expanded(z, ord, period)
  p <- length(e$ar)
  q <- length(e$ma)
  Jar <- matrix(0, p, P)
  Jma <- matrix(0, q, P)
  for (k in seq_len(nz)) {
    h <- 1e-6 * max(1, abs(z[k]))
    zp <- z; zp[k] <- zp[k] + h
    zm <- z; zm[k] <- zm[k] - h
    ep <- free_to_expanded(zp, ord, period)
    em <- free_to_expanded(zm, ord, period)
    if (p > 0L) Jar[, before + k] <- (ep$ar - em$ar) / (2 * h)
    if (q > 0L) Jma[, before + k] <- (ep$ma - em$ma) / (2 * h)
  }
  r <- arma_loglik_grad_cpp(as.numeric(x), if (is.null(dx)) matrix(0, 0, P) else dx, e$ar, Jar, e$ma, Jma, 1e-12)
  if (!is.finite(r$value) || any(!is.finite(r$gradient))) return(NULL)
  -r$gradient
}

# Gradient for optim(): exact where available, otherwise central differences of armacop_negll.
armacop_negll_gr <- function(z, ord, period, x) {
  g <- armacop_negll_grad(z, ord, period, x)
  if (!is.null(g)) return(g)
  vapply(seq_along(z), function(k) {
    h <- 1e-6
    zp <- z; zp[k] <- zp[k] + h
    zm <- z; zm[k] <- zm[k] - h
    (armacop_negll(zp, ord, period, x) - armacop_negll(zm, ord, period, x)) / (2 * h)
  }, numeric(1))
}

# Shrink AR-type coefficients until all partial autocorrelations are inside
# (-bound, bound), so a starting value is comfortably admissible.
arma_shrink <- function(ar, bound = 0.98) {
  if (length(ar) == 0L) return(ar)
  ok <- function(a) {
    pc <- ar_to_pacf(a)
    !anyNA(pc) && all(abs(pc) < bound)
  }
  for (i in 1:500) {
    if (ok(ar)) return(ar)
    ar <- ar * 0.95^seq_along(ar)
  }
  numeric(length(ar))
}

# Hannan-Rissanen starting values for the normal scores x: a long AR fit gives
# innovation estimates, then x is regressed on its own lags and the lagged
# innovations at the lags of the model (for seasonal models, the AR and MA
# lags and the seasonal lags; cross-terms of the multiplicative structure are
# ignored, as the optimiser refines the start). Returns a coefficient list.
arma_hr_start <- function(x, ord, period) {
  x[is.na(x)] <- 0
  n <- length(x)
  s <- if (is.na(period)) 0 else period
  arlags <- c(seq_len(ord[1L]), s * seq_len(ord[3L]))
  malags <- c(seq_len(ord[2L]), s * seq_len(ord[4L]))
  maxlag <- max(c(arlags, malags, 0L))
  k <- min(max(ceiling(2 * log(n)), 2 * s, maxlag + 1L, 2L), floor(n / 4))
  nreg <- length(arlags) + length(malags)
  first <- k + maxlag + 1L
  if (n - first + 1L < nreg + 10L) {
    stop("too few observations to compute starting values; supply 'init'")
  }
  e <- stats::ar.yw(x, aic = FALSE, order.max = k, demean = FALSE)$resid
  e[is.na(e)] <- 0
  t <- first:n
  X <- cbind(
    vapply(arlags, function(l) x[t - l], numeric(length(t))),
    vapply(malags, function(l) e[t - l], numeric(length(t)))
  )
  b <- stats::lm.fit(X, x[t])$coefficients
  b[is.na(b)] <- 0
  b <- unname(b)
  nar <- length(arlags)
  pars <- list(
    ar = b[seq_len(ord[1L])],
    sar = b[ord[1L] + seq_len(ord[3L])],
    ma = b[nar + seq_len(ord[2L])],
    sma = b[nar + ord[2L] + seq_len(ord[4L])]
  )
  pars$ar <- arma_shrink(pars$ar)
  pars$sar <- arma_shrink(pars$sar)
  pars$ma <- -arma_shrink(-pars$ma)
  pars$sma <- -arma_shrink(-pars$sma)
  pars[c("ar", "ma", "sar", "sma")]
}

# Starting model for tsfit: 'init' if given, else the coefficients of 'object'
# if it is a known model, else Hannan-Rissanen values.
armacop_start <- function(object, init, x) {
  if (!is.null(init)) {
    if (is(init, "armacop")) {
      if (!identical(init@order, object@order) || !identical(init@period, object@period)) {
        stop("'init' must have the same orders and period as 'object'")
      }
      if (!arma_is_specified(init)) stop("'init' has unspecified coefficients")
      return(init)
    }
    if (is.numeric(init) && length(init) == sum(object@order)) {
      pars <- split_blocks(as.numeric(init), object@order)
      return(armacop(ar = pars$ar, ma = pars$ma, sar = pars$sar, sma = pars$sma,
                     period = object@period))
    }
    stop("'init' must be an armacop model or a numeric vector of coefficients")
  }
  if (arma_is_specified(object)) return(object)
  pars <- arma_hr_start(x, object@order, object@period)
  armacop(ar = pars$ar, ma = pars$ma, sar = pars$sar, sma = pars$sma,
          period = object@period)
}

# Natural coefficients (ar, ma, sar, sma order) from the unconstrained vector.
free_to_coef <- function(z, ord) {
  b <- split_blocks(z, ord)
  c(free_to_ar(b$ar), free_to_ma(b$ma), free_to_ar(b$sar), free_to_ma(b$sma))
}

# Central-difference Hessian of f at z.
num_hessian <- function(f, z, h = 1e-3) {
  k <- length(z)
  H <- matrix(0, k, k)
  f0 <- f(z)
  e <- function(i) { v <- numeric(k); v[i] <- h; v }
  for (i in seq_len(k)) {
    H[i, i] <- (f(z + e(i)) - 2 * f0 + f(z - e(i))) / h^2
    for (j in seq_len(i - 1L)) {
      H[i, j] <- H[j, i] <- (f(z + e(i) + e(j)) - f(z + e(i) - e(j)) -
                               f(z - e(i) + e(j)) + f(z - e(i) - e(j))) / (4 * h^2)
    }
  }
  H
}

# Central-difference Jacobian of the vector function g at z.
num_jacobian <- function(g, z, h = 1e-6) {
  k <- length(z)
  J <- matrix(0, length(g(z)), k)
  for (i in seq_len(k)) {
    d <- numeric(k)
    d[i] <- h
    J[, i] <- (g(z + d) - g(z - d)) / (2 * h)
  }
  J
}

# Maximum likelihood on normal scores of copula data u. Returns the fitted
# armacop model, the starting model, the optimiser output and the optimum on
# the unconstrained scale.
armacop_mle <- function(u, object, init, control) {
  x <- stats::qnorm(u)
  ord <- object@order
  period <- object@period
  start <- armacop_start(object, init, x)
  z0 <- pars_to_free(start@pars)
  ctrl <- modifyList(
    list(maxit = 1000, reltol = 1e-12, ndeps = rep(1e-6, length(z0))), control
  )
  opt <- stats::optim(z0, armacop_negll, gr = armacop_negll_gr, ord = ord, period = period, x = x,
                      method = "BFGS", control = ctrl)
  if (opt$convergence != 0L) {
    alt <- stats::optim(opt$par, armacop_negll, ord = ord, period = period,
                        x = x, method = "Nelder-Mead",
                        control = list(maxit = 5000, reltol = 1e-12))
    if (alt$value < opt$value) opt <- alt
  }
  b <- split_blocks(opt$par, ord)
  model <- armacop(
    ar = free_to_ar(b$ar), ma = free_to_ma(b$ma),
    sar = free_to_ar(b$sar), sma = free_to_ma(b$sma), period = period
  )
  list(model = model, start = start, opt = opt, x = x)
}

# Hessian standard errors: inverse observed information on the unconstrained
# scale, carried to the coefficients by the delta method. Returns a covariance
# matrix, or NULL (with a warning) if the information matrix is not positive
# definite.
armacop_hessian_vcov <- function(fit, ord, period) {
  z <- fit$opt$par
  f <- function(v) armacop_negll(v, ord, period, fit$x)
  H <- num_hessian(f, z)
  V <- tryCatch(solve(H), error = function(e) NULL)
  if (is.null(V) || any(!is.finite(V)) || any(diag(V) <= 0) ||
      any(eigen(H, symmetric = TRUE, only.values = TRUE)$values <= 0)) {
    warning("the information matrix is not positive definite; standard errors not available")
    return(NULL)
  }
  J <- num_jacobian(function(v) free_to_coef(v, ord), z)
  V <- J %*% V %*% t(J)
  (V + t(V)) / 2
}

#' @describeIn tsfit Maximum likelihood fit of an [armacop-class] model.
#'
#' @param pseudo logical: `TRUE` (default) if `U` are pseudo-observations
#'   (ranks), `FALSE` if they are probability integral transforms from fitted
#'   marginal models. It matters for the bootstrap, which re-ranks simulated
#'   data when `pseudo = TRUE`.
#' @param se how to compute standard errors: `FALSE` (default, none),
#'   `"hessian"` or `"bootstrap"` (`TRUE` means `"bootstrap"`).
#' @param B number of bootstrap replicates (for `se = "bootstrap"`).
#' @param init optional starting values: an [armacop-class] model with the
#'   same orders and period as `object`, or a numeric vector of coefficients
#'   in the order `ar`, `ma`, `sar`, `sma`. If `NULL`, the coefficients of
#'   `object` are used if it is a known model, and otherwise Hannan-Rissanen
#'   starting values are computed from the data.
#' @param control list of control parameters for [stats::optim()].
#'
#' @details The model is a zero-mean, unit-variance Gaussian ARMA/SARMA for
#'   the normal scores `qnorm(U)`, so only the ARMA coefficients are
#'   estimated. Estimation is by BFGS on an unconstrained scale (hyperbolic
#'   tangent of the partial autocorrelations), on which every parameter value
#'   is stationary and invertible. If BFGS does not converge, Nelder-Mead is
#'   tried from the BFGS value and the better result is kept.
#'
#'   `se = "hessian"` gives the inverse observed information, carried to the
#'   coefficients by the delta method. It is quick but only approximate: for
#'   pseudo-observations it ignores the variability of the ranks, and for PIT
#'   data it ignores estimation of the margins. `se = "bootstrap"` is a
#'   parametric bootstrap: `B` series are simulated from the fitted model
#'   (exact stationary start, same missing-value pattern), ranked if
#'   `pseudo = TRUE`, and refitted with the same starting rules as the
#'   original fit (from `init` or `object` if given, otherwise from
#'   Hannan-Rissanen values). The covariance is that of the replicates; a
#'   replicate whose refit fails is recorded as `NA` and dropped.
#'
#' @return An object of class [armacopfit-class].
#' @export
#'
#' @examples
#' set.seed(1)
#' x <- as.numeric(arima.sim(list(ar = 0.5, ma = 0.4), n = 500, sd = 0.7))
#' U <- (rank(x) - 0.5) / length(x)
#' tsfit(U, armacop(order = c(1, 1)), se = "hessian")
setMethod("tsfit", signature(U = "ANY", object = "armacop"),
  function(U, object, pseudo = TRUE, se = FALSE, B = 200, init = NULL,
           control = list(), ...) {
    if (isTRUE(se)) se <- "bootstrap"
    if (!identical(se, FALSE) && !(is.character(se) && length(se) == 1L &&
                                   se %in% c("hessian", "bootstrap"))) {
      stop("'se' must be FALSE, \"hessian\" or \"bootstrap\"")
    }
    if (identical(se, "bootstrap") &&
        (!is.numeric(B) || length(B) != 1L || B < 10 || B != round(B))) {
      stop("'B' must be a whole number of at least 10")
    }
    u <- check_copula_data(U)
    ord <- object@order
    period <- object@period
    npar <- as.integer(sum(ord))
    fit <- armacop_mle(u, object, init, control)
    model <- fit$model
    est <- coef(model)
    cn <- names(est)
    V <- matrix(NA_real_, npar, npar, dimnames = list(cn, cn))
    boot <- matrix(NA_real_, 0L, npar, dimnames = list(NULL, cn))
    se_method <- "none"
    if (identical(se, "hessian")) {
      Vh <- armacop_hessian_vcov(fit, ord, period)
      if (!is.null(Vh)) {
        V[] <- Vh
        se_method <- "hessian"
      }
    } else if (identical(se, "bootstrap")) {
      n <- length(u)
      missing <- is.na(u)
      boot <- matrix(NA_real_, B, npar, dimnames = list(NULL, cn))
      for (b in seq_len(B)) {
        ub <- tssim(model, n)
        ub[missing] <- NA
        if (pseudo) {
          ub <- rank(ub, na.last = "keep", ties.method = "random") / (sum(!missing) + 1)
        }
        fb <- tryCatch(armacop_mle(ub, object, init, control)$model,
                       error = function(e) NULL)
        if (!is.null(fb)) boot[b, ] <- coef(fb)
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
    new("armacopfit",
      model = model, data = u, pseudo = pseudo, estimate = est,
      se = se_vec, vcov = V, se_method = se_method, boot = boot,
      loglik = -fit$opt$value - sum(stats::dnorm(fit$x, log = TRUE), na.rm = TRUE),
      nobs = as.integer(sum(!is.na(u))), npar = npar,
      convergence = as.integer(fit$opt$convergence),
      message = if (is.null(fit$opt$message)) "" else fit$opt$message,
      start = fit$start,
      template = if (inherits(U, c("ts", "zoo", "xts"))) U else NULL
    )
  }
)

#' @rdname tscop-methods
#' @export
setMethod("coef", "armacopfit", function(object) object@estimate)

#' @rdname tscop-methods
#' @export
setMethod("logLik", "armacopfit", function(object) {
  structure(object@loglik, df = object@npar, nobs = object@nobs, class = "logLik")
})

#' @rdname tscop-methods
#' @export
setMethod("vcov", "armacopfit", function(object, ...) object@vcov)

#' @rdname tscop-methods
#' @export
setMethod("show", "armacopfit", function(object) {
  cat("object class: armacopfit\n")
  cat("name: ", object@model@name, "\n", sep = "")
  if (!is.na(object@model@period)) cat("period: ", object@model@period, "\n", sep = "")
  cat("estimates:\n")
  if (object@se_method == "none") {
    print(object@estimate)
  } else {
    tab <- cbind(estimate = object@estimate, se = object@se)
    print(tab)
    cat("standard errors: ", object@se_method, "\n", sep = "")
  }
  cat("loglik: ", format(object@loglik, digits = 7), "  npar: ", object@npar,
      "  nobs: ", object@nobs, "\n", sep = "")
  if (object@convergence != 0L) {
    cat("warning: optimiser did not converge (code ", object@convergence, ")\n", sep = "")
  }
  invisible(object)
})

#' @describeIn tssim Simulation from the fitted model; `n` defaults to the
#'   length of the data.
#' @export
setMethod("tssim", "armacopfit", function(object, n, ...) {
  if (missing(n)) n <- length(object@data)
  tssim(object@model, n)
})

#' @describeIn ktau Kendall's tau of the pair copulas of the fitted model.
#' @export
setMethod("ktau", "armacopfit", function(object, lagmax, ...) {
  if (missing(lagmax)) lagmax <- 20
  ktau(object@model, lagmax)
})

#' @rdname resid-method
#' @export
setMethod("resid", "armacopfit", function(object, data, scale = "uniform", ...) {
  if (missing(data)) {
    data <- if (is.null(object@template)) object@data else object@template
  }
  armacop_resid(object@model, data, scale)
})

#' @rdname resid-method
#' @export
setMethod("residuals", "armacopfit", function(object, data, scale = "uniform", ...) {
  if (missing(data)) {
    data <- if (is.null(object@template)) object@data else object@template
  }
  armacop_resid(object@model, data, scale)
})

#' @rdname predict-method
#' @export
setMethod("predict", "armacopfit", function(object, data, h = 1, type = "moments",
                                            x = NULL, ...) {
  if (missing(data)) data <- object@data
  armacop_predict(object@model, data, h, type, x)
})
