#' ARMA copula processes
#'
#' The Gaussian copula of a stationary, invertible (seasonal) ARMA process with
#' mean zero and variance one. If \eqn{u_t} are the copula data, the normal
#' scores \eqn{x_t = \Phi^{-1}(u_t)} follow the ARMA process
#' \deqn{\Phi(B)\,\Phi_s(B^s)\, x_t = \Theta(B)\,\Theta_s(B^s)\, \epsilon_t,}
#' where the AR polynomial is \eqn{1 - \sum_i \mathrm{ar}_i B^i}, the MA
#' polynomial is \eqn{1 + \sum_j \mathrm{ma}_j B^j}, and the seasonal
#' polynomials are defined in the same way with `sar` and `sma` at multiples of
#' the period \eqn{s}. This is the sign convention of [stats::arima()].
#'
#' The innovation variance is not a parameter: it is determined by the
#' coefficients so that \eqn{\mathrm{Var}(x_t) = 1}. The free parameters are
#' exactly the ARMA coefficients.
#'
#' An object is either a *known model* (all coefficients given) or a
#' *structure to be fitted*, which records the orders and holds `NA`
#' coefficients.
#'
#' @slot name a short description, e.g. `"ARMA(1,1)"` or
#'   `"SARMA(1,0)(0,1)[12]"`.
#' @slot order integer vector `c(p, q, P, Q)` with the numbers of AR, MA,
#'   seasonal AR and seasonal MA coefficients.
#' @slot period the seasonal period (`NA` if there are no seasonal terms).
#' @slot pars a list with numeric components `ar`, `ma`, `sar` and `sma`
#'   (each possibly of length zero); `NA` entries mark a structure that has not
#'   been fitted.
#'
#' @exportClass armacop
setClass("armacop", contains = "gausscop", slots = c(
  name = "character",
  order = "integer",
  period = "numeric",
  pars = "list"
))

setValidity("armacop", function(object) {
  ord <- object@order
  if (length(ord) != 4L || anyNA(ord) || any(ord < 0L)) {
    return("'order' must be four non-negative integers (p, q, P, Q)")
  }
  pars <- object@pars
  if (!identical(names(pars), c("ar", "ma", "sar", "sma"))) {
    return("'pars' must be a list with components ar, ma, sar and sma")
  }
  lens <- vapply(pars, length, integer(1))
  if (any(lens != ord)) {
    return("the lengths of the coefficient vectors do not match 'order'")
  }
  if (!all(vapply(pars, is.numeric, logical(1)))) {
    return("coefficients must be numeric")
  }
  if (sum(ord) == 0L) {
    return("specify at least one ARMA coefficient")
  }
  period <- object@period
  if (length(period) != 1L) {
    return("'period' must be a single number")
  }
  if (ord[3L] + ord[4L] > 0L) {
    if (is.na(period) || period < 2 || period != round(period)) {
      return("'period' must be a whole number of at least 2 when seasonal terms are present")
    }
  } else if (!is.na(period)) {
    return("'period' must be NA when there are no seasonal terms")
  }
  if (!anyNA(unlist(pars))) {
    if (!all(is.finite(unlist(pars)))) {
      return("coefficients must be finite")
    }
    if (!arma_is_stationary(pars$ar)) return("non-stationary AR part")
    if (!arma_is_invertible(pars$ma)) return("non-invertible MA part")
    if (!arma_is_stationary(pars$sar)) return("non-stationary seasonal AR part")
    if (!arma_is_invertible(pars$sma)) return("non-invertible seasonal MA part")
  }
  TRUE
})

#' Constructor for ARMA copula processes
#'
#' Creates an [armacop-class] object, either as a known model by giving
#' the coefficients, or as a structure to be fitted by giving the orders.
#'
#' @param ar,ma numeric vectors of AR and MA coefficients (known model).
#' @param sar,sma numeric vectors of seasonal AR and MA coefficients (known
#'   model); these require `period`.
#' @param order integer vector `c(p, q)` giving the numbers of AR and MA
#'   coefficients of a structure to be fitted.
#' @param seasonal integer vector `c(P, Q)` giving the numbers of seasonal AR
#'   and MA coefficients of a structure to be fitted; requires `period` unless
#'   it is `c(0, 0)`.
#' @param period the seasonal period, e.g. 12 for monthly and 4 for quarterly
#'   data. It is ignored (stored as `NA`) when there are no seasonal terms.
#'
#' @details Give either coefficients (`ar`, `ma`, `sar`, `sma`) or orders
#'   (`order`, `seasonal`), not both. A structure has `NA` coefficients; it
#'   describes a model to be passed to a fitting function, which supplies
#'   starting values. The ARMA(0,0) model, independence, is not allowed.
#'
#' @return An object of class [armacop-class].
#' @export
#'
#' @examples
#' armacop(ar = 0.5, ma = 0.4)
#' armacop(ar = 0.5, sar = 0.3, sma = 0.2, period = 12)
#' armacop(order = c(1, 1))
#' armacop(order = c(1, 1), seasonal = c(0, 1), period = 12)
armacop <- function(ar = numeric(0), ma = numeric(0),
                       sar = numeric(0), sma = numeric(0),
                       order = NULL, seasonal = NULL, period = NA) {
  coefs <- list(ar = ar, ma = ma, sar = sar, sma = sma)
  if (any(vapply(coefs, is.null, logical(1)))) {
    stop("coefficients cannot be NULL; omit them or use numeric(0)")
  }
  structure_mode <- !is.null(order) || !is.null(seasonal)
  if (structure_mode) {
    if (any(lengths(coefs) > 0L)) {
      stop("give either coefficients or 'order' and 'seasonal', not both")
    }
    if (is.null(order)) stop("'order' is required for a model structure")
    if (is.null(seasonal)) seasonal <- c(0, 0)
    for (nm in c("order", "seasonal")) {
      v <- get(nm)
      if (!is.numeric(v) || length(v) != 2L || anyNA(v) || any(v < 0) ||
          any(v != round(v))) {
        stop(sprintf("'%s' must be two non-negative whole numbers", nm))
      }
    }
    ord <- as.integer(c(order, seasonal))
    coefs <- lapply(ord, function(k) rep(NA_real_, k))
    names(coefs) <- c("ar", "ma", "sar", "sma")
  } else {
    for (nm in names(coefs)) {
      v <- coefs[[nm]]
      if (!is.numeric(v) || anyNA(v)) {
        stop(sprintf("'%s' must be a numeric vector without missing values", nm))
      }
      coefs[[nm]] <- as.numeric(v)
    }
    ord <- as.integer(lengths(coefs))
  }
  seasonal_terms <- ord[3L] + ord[4L] > 0L
  if (seasonal_terms) {
    if (length(period) != 1L || is.na(period)) {
      stop("'period' is required when there are seasonal terms")
    }
    period <- as.numeric(period)
  } else {
    period <- NA_real_
  }
  for (nm in names(coefs)) {
    if (length(coefs[[nm]]) > 0L) {
      names(coefs[[nm]]) <- paste0(nm, seq_along(coefs[[nm]]))
    }
  }
  name <- if (seasonal_terms) {
    sprintf("SARMA(%d,%d)(%d,%d)[%g]", ord[1L], ord[2L], ord[3L], ord[4L], period)
  } else {
    sprintf("ARMA(%d,%d)", ord[1L], ord[2L])
  }
  new("armacop", name = name, order = ord, period = period, pars = coefs)
}

# TRUE if all coefficients are specified (a known model, or a fitted one).
arma_is_specified <- function(object) {
  !anyNA(unlist(object@pars))
}

# Expanded plain ARMA coefficients of a known model.
arma_coef_expanded <- function(object) {
  if (!arma_is_specified(object)) {
    stop("the model has unspecified coefficients; fit it to data first")
  }
  p <- lapply(object@pars, unname)
  arma_expand(p$ar, p$ma, p$sar, p$sma, object@period)
}

#' @describeIn armacop-class Coefficients, in the order `ar`, `ma`, `sar`,
#'   `sma`; `NA` for a structure that has not been fitted.
#' @param object an object of class [armacop-class].
#' @export
setMethod("coef", "armacop", function(object) {
  out <- unlist(unname(object@pars))
  if (is.null(out)) numeric(0) else out
})

#' @describeIn armacop-class Show method.
#' @export
setMethod("show", "armacop", function(object) {
  cat("object class: armacop\n")
  cat("name: ", object@name, "\n", sep = "")
  if (!is.na(object@period)) cat("period: ", object@period, "\n", sep = "")
  if (arma_is_specified(object)) {
    cat("parameters:\n")
  } else {
    cat("parameters (not yet estimated):\n")
  }
  print(coef(object))
  invisible(object)
})

#' @describeIn tssim Exact simulation from a known [armacop-class] model
#'   (`n` defaults to 1000). The
#'   first state is drawn from the stationary distribution, so no burn-in is
#'   needed. Returns `pnorm` of the simulated zero-mean, unit-variance ARMA
#'   process.
#' @export
setMethod("tssim", "armacop", function(object, n, ...) {
  if (missing(n)) n <- 1000
  if (!is.numeric(n) || length(n) != 1L || n < 1 || n != round(n)) {
    stop("'n' must be a positive whole number")
  }
  e <- arma_coef_expanded(object)
  stats::pnorm(arma_simulate(n, e$ar, e$ma))
})

# Partial autocorrelations at lags 1, ..., lag.max of the ARMA process, by the
# Durbin-Levinson recursion on its autocorrelations.
arma_pacf <- function(ar, ma, lag.max) {
  if (lag.max < 1L) return(numeric(0))
  g <- arma_acvf(ar, ma, lag.max)
  rho <- g[-1L] / g[1L]
  pacf <- numeric(lag.max)
  phi <- numeric(0)
  for (k in seq_len(lag.max)) {
    if (k == 1L) {
      pk <- rho[1L]
    } else {
      pk <- (rho[k] - sum(phi * rho[(k - 1L):1L])) / (1 - sum(phi * rho[1:(k - 1L)]))
    }
    pacf[k] <- pk
    phi <- c(phi - pk * rev(phi), pk)
  }
  pacf
}

#' @describeIn ktau Kendall's tau of the Gaussian pair copulas of an
#'   [armacop-class] model (`lagmax` defaults to 20).
#' @export
setMethod("ktau", "armacop", function(object, lagmax, ...) {
  if (missing(lagmax)) lagmax <- 20
  if (!is.numeric(lagmax) || length(lagmax) != 1L || lagmax < 1 ||
      lagmax != round(lagmax)) {
    stop("'lagmax' must be a positive whole number")
  }
  e <- arma_coef_expanded(object)
  (2 / pi) * asin(arma_pacf(e$ar, e$ma, lagmax))
})

# Restore ts/zoo/xts attributes of a template on a result of the same length.
restore_attrs <- function(out, template) {
  if (!is.null(template) && inherits(template, c("ts", "zoo", "xts")) &&
      length(out) == length(template)) {
    attributes(out) <- attributes(template)
  }
  out
}

# Standardised one-step prediction errors of the normal scores.
armacop_resid <- function(model, data, scale) {
  scale <- match.arg(scale, c("uniform", "normal"))
  u <- check_copula_data(data, min_n = 1L)
  e <- arma_coef_expanded(model)
  f <- arma_filter(stats::qnorm(u), e$ar, e$ma)
  r <- f$v / sqrt(f$F)
  out <- if (scale == "uniform") stats::pnorm(r) else r
  restore_attrs(out, data)
}

#' Residuals of ARMA copula models
#'
#' The residuals are the standardised one-step prediction errors of the normal
#' scores \eqn{\Phi^{-1}(u_t)}. If the model is correct they are independent
#' standard normal, so on the default `"uniform"` scale (their `pnorm`) they
#' are iid uniform. Missing observations give `NA`. If the data are a `ts`,
#' `zoo` or `xts` object the time index is restored.
#'
#' @param object an [armacop-class] model with known coefficients, or an
#'   [armacopfit-class] object.
#' @param data copula data in (0, 1) for a known model; for a fitted object it
#'   defaults to the data to which the model was fitted.
#' @param scale `"uniform"` (default) or `"normal"`.
#' @param ... not used.
#'
#' @return A numeric vector the length of the data.
#' @name resid-armacop
#' @aliases resid,armacop-method residuals,armacop-method
NULL

#' @rdname resid-armacop
#' @export
setMethod("resid", "armacop", function(object, data, scale = "uniform", ...) {
  armacop_resid(object, data, scale)
})

#' @rdname resid-armacop
#' @export
setMethod("residuals", "armacop", function(object, data, scale = "uniform", ...) {
  armacop_resid(object, data, scale)
})

#' Prediction for ARMA copula models
#'
#' Conditional distribution of the copula variable `h` steps ahead given the
#' past data. On the normal scale, \eqn{x_{n+h} = \Phi^{-1}(u_{n+h})} given
#' the past is normal with the mean and standard deviation returned by
#' `type = "moments"`. The other types transform this to the uniform scale.
#'
#' @param object an [armacop-class] model with known coefficients, or an
#'   [armacopfit-class] object.
#' @param data past copula data in (0, 1); for a fitted object it defaults to
#'   the data to which the model was fitted.
#' @param h horizon (number of steps ahead). For `type = "moments"` a vector
#'   of horizons is allowed; for the other types a single value.
#' @param type `"moments"` (default) for a data frame with columns `h`,
#'   `mean` and `sd`; `"df"` for the conditional distribution function of
#'   \eqn{u_{n+h}} at `x`; `"qf"` for its quantile function at probabilities
#'   `x`; `"dens"` for its density at `x`.
#' @param x points in (0, 1) at which to evaluate `"df"`, `"qf"` or `"dens"`.
#' @param ... not used.
#'
#' @return A data frame (`type = "moments"`) or a numeric vector the length of
#'   `x`.
#' @name predict-armacop
#' @aliases predict,armacop-method
NULL

# Validate the arguments of a Gaussian-process predict method; returns the
# matched type.
gauss_predict_check <- function(h, type, x) {
  type <- match.arg(type, c("moments", "df", "qf", "dens"))
  if (!is.numeric(h) || length(h) < 1L || anyNA(h) || any(h < 1) || any(h != round(h))) {
    stop("'h' must be positive whole numbers")
  }
  if (type != "moments") {
    if (length(h) != 1L) stop("a single horizon 'h' is required unless type = \"moments\"")
    if (is.null(x) || !is.numeric(x) || anyNA(x) || any(x <= 0 | x >= 1)) {
      stop("'x' must be numeric values strictly inside (0, 1)")
    }
  }
  type
}

# Predictive output from the conditional means mu[k] and variances v[k] of the
# normal score k steps ahead.
gauss_predict_result <- function(mu, v, h, type, x) {
  m <- mu[h]
  s <- sqrt(v[h])
  switch(type,
    moments = data.frame(h = h, mean = m, sd = s),
    df = stats::pnorm((stats::qnorm(x) - m) / s),
    qf = stats::pnorm(m + s * stats::qnorm(x)),
    dens = stats::dnorm(stats::qnorm(x), m, s) / stats::dnorm(stats::qnorm(x))
  )
}

armacop_predict <- function(model, data, h, type, x) {
  type <- gauss_predict_check(h, type, x)
  u <- check_copula_data(data, min_n = 1L)
  e <- arma_coef_expanded(model)
  ss <- arma_state_space(e$ar, e$ma)
  f <- arma_filter(stats::qnorm(u), e$ar, e$ma)
  a <- f$a
  P <- f$P
  hmax <- max(h)
  mu <- numeric(hmax)
  v <- numeric(hmax)
  for (k in seq_len(hmax)) {
    mu[k] <- a[1L]
    v[k] <- P[1L, 1L]
    a <- ss$phi * a[1L] + c(a[-1L], 0)
    P <- ss$Tm %*% P %*% t(ss$Tm) + ss$RR
  }
  gauss_predict_result(mu, v, h, type, x)
}

#' @rdname predict-armacop
#' @export
setMethod("predict", "armacop", function(object, data, h = 1, type = "moments",
                                         x = NULL, ...) {
  armacop_predict(object, data, h, type, x)
})
