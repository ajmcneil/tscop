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
