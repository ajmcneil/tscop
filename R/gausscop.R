#' Gaussian copula processes
#'
#' Virtual class of the Gaussian copula processes whose partial
#' autocorrelations can serve as the Kendall partial autocorrelation function
#' (kpacf) of a [dvinecop-class] model: a stationary Gaussian process with
#' standard normal margins is determined by its partial autocorrelations
#' \eqn{\alpha_k}, and the kpacf is \eqn{\tau_k = (2/\pi)\arcsin\alpha_k}.
#' Subclasses are [armacop-class] (ARMA and seasonal ARMA autocorrelations)
#' and `arfactcop` (a sum of independent AR(1) factors plus white noise,
#' with non-negative partial autocorrelations).
#'
#' @exportClass gausscop
setClass("gausscop", representation("VIRTUAL"))

setClassUnion("gausscopOrNULL", c("gausscop", "NULL"))

# ---- internal interface of a Gaussian copula process ---------------------
#
# Every subclass implements these generics; dvinecop, tsfit for dvinecop and
# autodvinecop use only this interface.
#
#   gc_is_specified(G)        TRUE if all coefficients are specified
#   gc_describe(G)            one-line description, e.g. "ARMA(1,1)"
#   gc_pacf(G, K)             partial autocorrelations at lags 1..K
#   gc_effective_lag(G, K0, tautol, maxlag)
#                             list(K, limit): effective lag of the kpacf
#   gc_structure(G)           a list identifying the model structure (orders)
#   gc_prefit(G, U, control)  Gaussian copula MLE of the structure/model G to
#                             copula data U: a fitted (known) model
#   gc_from_coef(G, v)        a model of the same structure with natural
#                             coefficient vector v
#   gc_free_ops(G)            closures for fitting on an unconstrained scale:
#                             list(nfree, zbound, to_free(G2), pacf(z, K),
#                             effective_lag(z, K0, tautol, maxlag), coef(z),
#                             build(z))

setGeneric("gc_is_specified", function(G) standardGeneric("gc_is_specified"))
setGeneric("gc_describe", function(G) standardGeneric("gc_describe"))
setGeneric("gc_pacf", function(G, K) standardGeneric("gc_pacf"))
setGeneric("gc_effective_lag", function(G, K0, tautol, maxlag) standardGeneric("gc_effective_lag"))
setGeneric("gc_structure", function(G) standardGeneric("gc_structure"))
setGeneric("gc_prefit", function(G, U, control) standardGeneric("gc_prefit"))
setGeneric("gc_from_coef", function(G, v) standardGeneric("gc_from_coef"))
setGeneric("gc_free_ops", function(G) standardGeneric("gc_free_ops"))

# Kendall's tau of the kpacf of G at lags 1..K.
gc_tau <- function(G, K) (2 / pi) * asin(gc_pacf(G, K))
