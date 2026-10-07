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
#'
#' @exportClass armacopfit
setClass("armacopfit", slots = c(
  model = "armacop", data = "numeric", pseudo = "logical",
  estimate = "numeric", se = "numeric", vcov = "matrix", se_method = "character",
  boot = "matrix", loglik = "numeric", nobs = "integer", npar = "integer",
  convergence = "integer", message = "character", start = "armacop"
))

# Validate copula data: numeric, in (0, 1), NA allowed.
check_copula_data <- function(U) {
  if (!is.numeric(U)) stop("'U' must be numeric")
  u <- as.numeric(U)
  ok <- !is.na(u)
  if (any(!is.finite(u[ok]))) stop("'U' must be finite")
  if (any(u[ok] <= 0 | u[ok] >= 1)) {
    stop("'U' must lie strictly inside (0, 1); see udp::boundaryadjust() for data on the boundary")
  }
  if (sum(ok) < 20L) stop("at least 20 non-missing observations are required")
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
  if (any(!is.finite(z)) || max(abs(z)) > 6) return(1e10)
  e <- free_to_expanded(z, ord, period)
  val <- tryCatch(-arma_filter(x, e$ar, e$ma)$loglik, error = function(err) NA_real_)
  if (!is.finite(val)) 1e10 else val
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

#' @describeIn tsfit Maximum likelihood fit of an [armacop-class] model.
#'
#' @param pseudo logical: `TRUE` (default) if `U` are pseudo-observations
#'   (ranks), `FALSE` if they are probability integral transforms from fitted
#'   marginal models. It matters for standard errors.
#' @param se how to compute standard errors. Only `FALSE` (none) is
#'   implemented so far.
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
#' @return An object of class [armacopfit-class].
#' @export
#'
#' @examples
#' set.seed(1)
#' x <- as.numeric(arima.sim(list(ar = 0.5, ma = 0.4), n = 500, sd = 0.7))
#' U <- (rank(x) - 0.5) / length(x)
#' tsfit(U, armacop(order = c(1, 1)))
setMethod("tsfit", signature(U = "ANY", object = "armacop"),
  function(U, object, pseudo = TRUE, se = FALSE, B = 200, init = NULL,
           control = list(), ...) {
    if (!identical(se, FALSE)) {
      stop("standard errors are not implemented yet; use se = FALSE")
    }
    u <- check_copula_data(U)
    x <- stats::qnorm(u)
    ord <- object@order
    period <- object@period
    start <- armacop_start(object, init, x)
    z0 <- pars_to_free(start@pars)
    ctrl <- utils::modifyList(
      list(maxit = 1000, reltol = 1e-12, ndeps = rep(1e-6, length(z0))), control
    )
    opt <- stats::optim(z0, armacop_negll, ord = ord, period = period, x = x,
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
    npar <- as.integer(sum(ord))
    new("armacopfit",
      model = model, data = u, pseudo = pseudo, estimate = coef(model),
      se = rep(NA_real_, npar), vcov = matrix(NA_real_, npar, npar),
      se_method = "none", boot = matrix(NA_real_, 0L, npar),
      loglik = -opt$value - sum(stats::dnorm(x, log = TRUE), na.rm = TRUE),
      nobs = as.integer(sum(!is.na(u))), npar = npar,
      convergence = as.integer(opt$convergence),
      message = if (is.null(opt$message)) "" else opt$message,
      start = start
    )
  }
)

#' @describeIn armacopfit-class Coefficient estimates.
#' @param object an object of class [armacopfit-class].
#' @export
setMethod("coef", "armacopfit", function(object) object@estimate)

#' @describeIn armacopfit-class Log-likelihood, with `df` and `nobs`
#'   attributes so that [stats::AIC()] and [stats::BIC()] work.
#' @export
setMethod("logLik", "armacopfit", function(object) {
  structure(object@loglik, df = object@npar, nobs = object@nobs, class = "logLik")
})

#' @describeIn armacopfit-class Show method.
#' @export
setMethod("show", "armacopfit", function(object) {
  cat("object class: armacopfit\n")
  cat("name: ", object@model@name, "\n", sep = "")
  if (!is.na(object@model@period)) cat("period: ", object@model@period, "\n", sep = "")
  cat("estimates:\n")
  print(object@estimate)
  cat("loglik: ", format(object@loglik, digits = 7), "  npar: ", object@npar,
      "  nobs: ", object@nobs, "\n", sep = "")
  if (object@convergence != 0L) {
    cat("warning: optimiser did not converge (code ", object@convergence, ")\n", sep = "")
  }
  invisible(object)
})
