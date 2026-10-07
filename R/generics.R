#' Fit a time series copula model
#'
#' Generic for maximum likelihood estimation of the time series copula models
#' in tscop. Methods exist for [armacop-class] objects.
#'
#' @param U a numeric vector or time series of copula data in (0, 1):
#'   pseudo-observations (ranks), or probability integral transforms from
#'   fitted marginal models (set `pseudo = FALSE`). `NA` values are allowed.
#'   Values equal to 0 or 1 are rejected; see [udp::boundaryadjust()].
#' @param object the model to fit: a model structure, or a known model whose
#'   coefficients are used as starting values.
#' @param ... further arguments passed to methods.
#'
#' @return A fitted model object.
#' @export
setGeneric("tsfit", function(U, object, ...) standardGeneric("tsfit"))

#' Simulate from a time series copula model
#'
#' Generic for simulation from the time series copula models in tscop.
#' Methods exist for [armacop-class] and [armacopfit-class] objects.
#'
#' @param object a model object.
#' @param n length of the simulated series.
#' @param ... further arguments passed to methods.
#'
#' @return A numeric vector of copula data in (0, 1) of length `n`.
#' @export
setGeneric("tssim", function(object, n, ...) standardGeneric("tssim"),
           signature = "object")

#' Kendall's tau of the pair copulas of a time series copula model
#'
#' Generic returning Kendall's tau of the pair copulas at lags
#' `1, ..., lagmax` in the D-vine representation of the process. For an
#' [armacop-class] model the pair copulas are Gaussian with parameters equal
#' to the partial autocorrelations of the ARMA process, so
#' \eqn{	au_k = (2/\pi)\arcsin(\alpha_k)} with \eqn{\alpha_k} the lag-`k`
#' partial autocorrelation.
#'
#' @param object a model object.
#' @param lagmax maximum lag.
#' @param ... further arguments passed to methods.
#'
#' @return A numeric vector of length `lagmax`.
#' @export
setGeneric("ktau", function(object, lagmax = 20, ...) standardGeneric("ktau"),
           signature = "object")
