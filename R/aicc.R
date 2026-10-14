# The corrected Akaike information criterion.

# AICc from the log-likelihood, the number of parameters and the number of observations; infinite
# when the correction is not defined (nobs <= npar + 1)
ic_aicc <- function(loglik, npar, nobs) {
  den <- nobs - npar - 1
  ifelse(den > 0, -2 * loglik + 2 * npar + 2 * npar * (npar + 1) / den, Inf)
}

#' Akaike corrected information criterion
#'
#' \eqn{\mathrm{AICc} = -2\ell + 2k + 2k(k + 1)/(n - k - 1)} for a model with
#' log-likelihood \eqn{\ell}, \eqn{k} estimated parameters and \eqn{n}
#' observations: the criterion AIC with the small-sample correction of Hurvich
#' and Tsai (1989). It works for every fitted model object of this package, like
#' [stats::AIC()] and [stats::BIC()] (which use the same `logLik` method), so that
#' models of all classes can be compared by AIC, BIC and AICc. It is `Inf` when
#' \eqn{n \le k + 1}. As for the other criteria, compare only models for the same
#' data and the same likelihood: a [tsmod-class] model with an empirical margin has a
#' pseudo-likelihood that is not comparable with the full likelihood of one with a parametric
#' margin.
#'
#' @param object a fitted model object with a `logLik` method.
#' @param ... further fitted model objects (for the same data): the result is then a
#'   data frame with a row per model, as for [stats::AIC()].
#'
#' @return A number, or for several objects a data frame with the columns `df` (the
#'   number of estimated parameters) and `AICc`.
#' @export
#' @examples
#' set.seed(1)
#' u <- tssim(armacop(ar = 0.6, ma = 0.3), 300)
#' f1 <- tsfit(u, armacop(order = c(1, 0)), pseudo = FALSE)
#' f2 <- tsfit(u, armacop(order = c(1, 1)), pseudo = FALSE)
#' AICc(f1)
#' AICc(f1, f2)
#' c(AIC = AIC(f2), BIC = BIC(f2), AICc = AICc(f2))
AICc <- function(object, ...) UseMethod("AICc")

#' @export
AICc.default <- function(object, ...) {
  one <- function(ll) {
    no <- attr(ll, "nobs")
    if (is.null(no)) stop("the log-likelihood has no 'nobs' attribute", call. = FALSE)
    k <- attr(ll, "df")
    c(df = k, AICc = ic_aicc(as.numeric(ll), k, no))
  }
  if (...length() == 0L) return(unname(one(logLik(object))["AICc"]))
  objs <- list(object, ...)
  vals <- t(vapply(objs, function(o) one(logLik(o)), numeric(2)))
  val <- data.frame(df = vals[, "df"], AICc = vals[, "AICc"])
  Call <- match.call()
  row.names(val) <- as.character(Call[-1L])
  val
}
