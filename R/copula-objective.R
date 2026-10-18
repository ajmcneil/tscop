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
#   grad       NULL, or function(z, u, du) returning the gradient of the negative
#              log-likelihood with respect to nm directions along which the data u
#              vary (the columns of the n x nm matrix du) followed by z; NULL when
#              the exact gradient is unavailable
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
    K = NA_integer_,
    grad = function(z, u, du) {
      nx <- ncol(du)
      x <- stats::qnorm(u)
      dx <- du / stats::dnorm(x)
      g <- armacop_negll_grad(z, ord, period, x, cbind(dx, matrix(0, nrow(dx), length(z))), before = nx)
      if (is.null(g)) return(NULL)
      g[seq_len(nx)] <- g[seq_len(nx)] - colSums(x * dx)    # + sum log phi(x_t)
      g
    }
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
    K = NA_integer_,
    grad = function(z, u, du) {
      nx <- ncol(du)
      x <- stats::qnorm(u)
      dx <- du / stats::dnorm(x)                      # dx/dtheta = du/dtheta / phi(x)
      g <- afc_negll_grad(z, m, noise, x, cbind(dx, matrix(0, nrow(dx), length(z))), before = nx)
      if (is.null(g)) return(NULL)
      # + sum log phi(x_t): derivative -sum x_t dx_t along the data directions
      g[seq_len(nx)] <- g[seq_len(nx)] - colSums(x * dx)
      g
    }
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
    K = K,
    grad = function(z, u, du) {
      nx <- ncol(du)
      if (nx + length(z) > 12L) return(NULL)
      g <- dv_grad_pieces(z, lay, K)
      if (is.null(g)) return(NULL)
      e <- g$e
      pad <- function(m) cbind(matrix(0, nrow(m), nx), m)
      r <- dv_loglik_grad_cpp(u, cbind(du, matrix(0, nrow(du), length(z))), e$fam, e$rot, e$par, pad(g$dpar), e$ivt, e$bfam, e$brot,
                              e$bpar, pad(g$dbpar), e$d1, e$d2, c(numeric(nx), g$dd1), c(numeric(nx), g$dd2),
                              e$nl, e$k1, e$k2, c(numeric(nx), g$dk1), c(numeric(nx), g$dk2))
      if (!is.finite(r$value) || any(!is.finite(r$gradient))) return(NULL)
      -r$gradient
    }
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
    K = K,
    grad = function(z, u, du) vt_negll_grad(z, L, K, u, du)
  )
})
