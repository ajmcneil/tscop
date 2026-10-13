# Internal copula-objective interface for joint maximum likelihood.
#
# cop_obj(model, u) describes the copula likelihood of a fully specified fitted
# copula process on the unconstrained scale, independently of how the class
# fits itself:
#
#   z0         the free-scale vector of the model's coefficients
#   negll(z,u) negative copula log-likelihood of PITs u at the free vector z
#              (1e10 when infeasible or on numerical failure)
#   build(z)   the copula model with free vector z
#   coef(z)    the natural coefficient vector (named) at z
#   K          effective lag held fixed in the objective (NA if none)
#
# A joint fit of a margin and a copula process uses only this interface, so a
# class gets joint maximum likelihood by implementing one method (and exact
# gradients can later replace the numerical ones behind the same interface).

setGeneric("cop_obj", function(model, u) standardGeneric("cop_obj"))

setMethod("cop_obj", "armacop", function(model, u) {
  ops <- gc_free_ops(model)
  ord <- model@order
  period <- model@period
  list(
    z0 = ops$to_free(model),
    negll = function(z, u) {
      x <- stats::qnorm(u)
      v <- armacop_negll(z, ord, period, x)
      if (v >= 1e10) 1e10 else v + sum(stats::dnorm(x, log = TRUE))
    },
    build = ops$build,
    coef = function(z) stats::setNames(ops$coef(z), names(coef(model))),
    K = NA_integer_
  )
})

setMethod("cop_obj", "arfactcop", function(model, u) {
  ops <- gc_free_ops(model)
  m <- model@m
  noise <- model@noise
  list(
    z0 = ops$to_free(model),
    negll = function(z, u) {
      x <- stats::qnorm(u)
      v <- afc_negll(z, m, noise, x)
      if (v >= 1e10) 1e10 else v + sum(stats::dnorm(x, log = TRUE))
    },
    build = ops$build,
    coef = function(z) stats::setNames(ops$coef(z), names(coef(model))),
    K = NA_integer_
  )
})

setMethod("cop_obj", "sdvinecop", function(model, u) {
  if (model@mode != "tied") stop("joint estimation needs a tied sdvinecop model", call. = FALSE)
  st <- dv_start(model, NULL, u, list())
  lay <- st$lay
  K <- model@K
  limit <- model@Klimit
  list(
    z0 = st$z0,
    negll = function(z, u) dv_negll(z, lay, K, u),
    build = function(z) dv_build_model(model, z, lay, K, limit),
    coef = function(z) stats::setNames(dv_free_to_coef(z, lay), names(coef(model))),
    K = K
  )
})

setMethod("cop_obj", "vtscop", function(model, u) {
  st <- vt_start(model, NULL, u, list())
  L <- st$L
  base <- model@base
  K <- if (L$kind == "afc") NA_integer_ else base@K
  limit <- if (L$kind == "afc") "none" else base@Klimit
  list(
    z0 = st$z0,
    negll = function(z, u) vt_negll(z, L, K, u),
    build = function(z) vt_build_model(model, z, L, K, limit),
    coef = function(z) stats::setNames(vt_free_to_coef(z, L), names(coef(model))),
    K = K
  )
})
