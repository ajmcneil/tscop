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
#' @exportClass armacopula
setClass("armacopula", slots = c(
  name = "character",
  order = "integer",
  period = "numeric",
  pars = "list"
))

setValidity("armacopula", function(object) {
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
#' Creates an [armacopula-class] object, either as a known model by giving
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
#' @return An object of class [armacopula-class].
#' @export
#'
#' @examples
#' armacopula(ar = 0.5, ma = 0.4)
#' armacopula(ar = 0.5, sar = 0.3, sma = 0.2, period = 12)
#' armacopula(order = c(1, 1))
#' armacopula(order = c(1, 1), seasonal = c(0, 1), period = 12)
armacopula <- function(ar = numeric(0), ma = numeric(0),
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
  new("armacopula", name = name, order = ord, period = period, pars = coefs)
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

#' @describeIn armacopula-class Coefficients, in the order `ar`, `ma`, `sar`,
#'   `sma`; `NA` for a structure that has not been fitted.
#' @param object an object of class [armacopula-class].
#' @export
setMethod("coef", "armacopula", function(object) {
  out <- unlist(unname(object@pars))
  if (is.null(out)) numeric(0) else out
})

#' @describeIn armacopula-class Show method.
#' @export
setMethod("show", "armacopula", function(object) {
  cat("object class: armacopula\n")
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
