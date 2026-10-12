#' Process-level v-transform copula processes
#'
#' A copula process obtained from a base copula process \eqn{(V_t)} by
#' stochastic inversion of a v-transform \eqn{\mathcal V}: \eqn{U_t =
#' \mathcal V^{\leftarrow}(V_t, W_t)} with independent uniform \eqn{W_t}
#' (McNeil, 2021; Bladt and McNeil, 2021). The copula density of
#' \eqn{(U_1, \ldots, U_n)} is \eqn{c_V(\mathcal V(u_1), \ldots, \mathcal
#' V(u_n))}, where \eqn{c_V} is the density of the base process, for **any**
#' v-transform, so the likelihood is the base likelihood evaluated at the
#' transformed data and no h-function of the v-transform is ever needed.
#' Simulation applies the stochastic inverse to a simulated base series.
#'
#' The base process has non-negative Kendall partial autocorrelations: it is an
#' [arfactcop-class] model (the Gaussian copula of independent AR(1) factors
#' plus noise, a generalisation of the VT-ARMA model of McNeil, 2021) or a
#' tied [sdvinecop-class] model whose kpacf is an `arfactcop` (pair copulas with
#' non-Gaussian dependence, without inverse-v-transformed pair copulas).
#'
#' The v-transform is any of the families of the udp package: `vsymmetric()`
#' (fixed fulcrum 0.5), `vlinear(delta)`, `v2p(delta, kappa)`,
#' `v2b(delta, kappa)`, `v3p(delta, kappa, xi)` and `v3b(delta, kappa, xi)`.
#' Its fulcrum `delta` appears in the coefficients under that name (`delta1`
#' and `delta2` in the pairwise inverse-v-transformed models).
#'
#' @slot name short description.
#' @slot base the base process: an [arfactcop-class] or a tied [sdvinecop-class]
#'   object (specified or a structure to be fitted).
#' @slot vt a v-transform object of the udp package; its parameter values are
#'   the known values or the starting values.
#' @slot vtfix names of v-transform parameters (`"delta"`, `"kappa"`,
#'   `"xi"`) held fixed in fitting.
#'
#' @references McNeil, A. J. (2021). Modelling volatile time series with
#'   v-transforms and copulas. *Risks*, 9(1), 14.
#'
#'   Bladt, M. and McNeil, A. J. (2021). Time series copula models using
#'   d-vines and v-transforms. *Econometrics and Statistics*.
#'
#' @exportClass vtscop
setClass("vtscop", slots = c(name = "character", base = "ANY", vt = "vtransform",
                             vtfix = "character"))

# Parameters (delta, kappa, xi) of a udp v-transform.
vt_params <- function(vt) vt@pars

# A udp v-transform from its family name and tscop-named parameters.
vt_build <- function(name, p) {
  switch(name,
    vsymmetric = udp::vsymmetric(),
    vlinear = udp::vlinear(p[["delta"]]),
    v2p = udp::v2p(p[["delta"]], p[["kappa"]]),
    v2b = udp::v2b(p[["delta"]], p[["kappa"]]),
    v3p = udp::v3p(p[["delta"]], p[["kappa"]], p[["xi"]]),
    v3b = udp::v3b(p[["delta"]], p[["kappa"]], p[["xi"]]),
    stop("unsupported v-transform '", name, "'")
  )
}

# Check the base process of a vtscop; returns NULL if acceptable, else a message.
vt_check_base <- function(base) {
  if (is(base, "arfactcop")) return(NULL)
  if (is(base, "sdvinecop")) {
    if (base@mode != "tied") return("a sdvinecop base must be a tied model")
    if (!is(base@kpacf, "arfactcop")) {
      return("the kpacf of a sdvinecop base must be an arfactcop (non-negative Kendall partial autocorrelations)")
    }
    if (dv_has_ivt(base)) return("a sdvinecop base cannot contain inverse-v-transformed pair copulas")
    return(NULL)
  }
  "the base must be an arfactcop or a tied sdvinecop with an arfactcop kpacf"
}

setValidity("vtscop", function(object) {
  msg <- vt_check_base(object@base)
  if (!is.null(msg)) return(msg)
  free <- names(vt_params(object@vt))
  if (!all(object@vtfix %in% free)) {
    return(paste0("'vtfix' must name parameters of the v-transform: ",
                  if (length(free)) paste(free, collapse = ", ") else "(none)"))
  }
  TRUE
})

#' Constructor for process-level v-transform copula processes
#'
#' @param base the base copula process: an [arfactcop-class] model, or a tied
#'   [sdvinecop-class] model whose kpacf is an `arfactcop` and which has no
#'   inverse-v-transformed pair copula; specified, or a structure to be fitted.
#' @param vt a v-transform from the udp package: `udp::vsymmetric()`,
#'   `udp::vlinear(delta)`, `udp::v2p(delta, kappa)`, `udp::v2b(delta, kappa)`,
#'   `udp::v3p(delta, kappa, xi)` or `udp::v3b(delta, kappa, xi)`. The parameter
#'   values are the known values, or starting values for fitting (`delta` is
#'   the fulcrum).
#' @param vtfix character vector of v-transform parameters (`"delta"`,
#'   `"kappa"`, `"xi"`) that are held fixed at their values in `vt` when the
#'   model is fitted.
#'
#' @return An object of class [vtscop-class].
#' @export
#'
#' @examples
#' # VT-ARMA-type model: one AR factor plus noise, asymmetric v-transform
#' vtscop(arfactcop(lambda = 0.9, w = 0.6), udp::v2p(delta = 0.45, kappa = 0.8))
#' # a d-vine base with a Gumbel copula at lag 1
#' vtscop(sdvinecop(arfactcop(lambda = 0.9, w = 0.6), pairs = list(pc("gumbel"))),
#'        udp::vlinear(0.5))
#' # a model structure to be fitted, with the fulcrum fixed
#' vtscop(arfactcop(order = 2), udp::v2p(0.5, 1), vtfix = "delta")
vtscop <- function(base, vt, vtfix = character(0)) {
  msg <- vt_check_base(base)
  if (!is.null(msg)) stop(msg)
  if (!is(vt, "vtransform")) stop("'vt' must be a v-transform of the udp package (vlinear(), v2p(), ...)")
  vtfix <- as.character(vtfix)
  bname <- if (is(base, "sdvinecop")) "d-vine" else base@name
  name <- paste0("VT(", vt@name, ", ", bname, ")")
  new("vtscop", name = name, base = base, vt = vt, vtfix = vtfix)
}

# Is the base specified?
vt_base_specified <- function(object) {
  b <- object@base
  if (is(b, "arfactcop")) gc_is_specified(b) else !is.na(b@K) && gc_is_specified(b@kpacf)
}

# The base process's copula log-likelihood of v-transformed data v in (0, 1) at
# a specified model.
vt_base_loglik <- function(base, v) {
  if (is(base, "arfactcop")) {
    p <- afc_check_known(base)
    afc_copula_loglik(v, p$lambda, p$w, p$w0)
  } else {
    dv_model_loglik(base, v)
  }
}

# Clamped v-transformed data.
vt_transform_data <- function(vt, u) {
  v <- udp::udptrans(vt, u)
  pmin(pmax(v, 1e-10), 1 - 1e-10)
}

# Copula log-likelihood of complete or (for Gaussian bases) incomplete data u
# under a specified vtscop.
vt_loglik <- function(object, u) {
  ok <- !is.na(u)
  v <- rep(NA_real_, length(u))
  v[ok] <- vt_transform_data(object@vt, u[ok])
  vt_base_loglik(object@base, v)
}

#' @describeIn vtscop-class Coefficients: those of the base process, then the
#'   free parameters of the v-transform (`delta`, `kappa`, `xi`).
#' @param object an object of class [vtscop-class].
#' @export
setMethod("coef", "vtscop", function(object) {
  p <- vt_params(object@vt)
  free <- setdiff(names(p), object@vtfix)
  c(coef(object@base), p[free])
})

#' @describeIn vtscop-class Show method.
#' @export
setMethod("show", "vtscop", function(object) {
  cat("object class: vtscop\n")
  cat("v-transform: ", object@vt@name,
      if (length(object@vt@pars) > 0L) paste0(" (", paste(names(vt_params(object@vt)), "=",
                                                        format(vt_params(object@vt), digits = 4), collapse = ", "), ")"),
      if (length(object@vtfix) > 0L) paste0("; fixed: ", paste(object@vtfix, collapse = ", ")), "\n", sep = "")
  cat("base process:\n")
  b <- utils::capture.output(show(object@base))
  cat(paste0("  ", b), sep = "\n")
  invisible(object)
})

# Simulation of n values, nrep replicates (columns): base simulation followed by
# the stochastic inverse of the v-transform.
vt_simulate <- function(object, n, nrep = 1L) {
  b <- object@base
  V <- if (is(b, "arfactcop")) {
    p <- afc_check_known(b)
    vapply(seq_len(nrep), function(i) stats::pnorm(afc_simulate(n, p$lambda, p$w, p$w0)), numeric(n))
  } else {
    pairs <- dv_pairs_or_stop(b)
    dv_simulate(pairs, matrix(stats::runif(n * nrep), n, nrep))
  }
  V <- matrix(V, n, nrep)
  U <- udp::udpsi(object@vt, as.numeric(V), stats::runif(n * nrep))
  matrix(pmin(pmax(U, 1e-12), 1 - 1e-12), n, nrep)
}

#' @describeIn tssim Simulation from a known [vtscop-class] model: the base
#'   process is simulated and the stochastic inverse of the v-transform is
#'   applied (`n` defaults to 1000).
#' @export
setMethod("tssim", "vtscop", function(object, n, ...) {
  if (missing(n)) n <- 1000
  if (!is.numeric(n) || length(n) != 1L || n < 1 || n != round(n)) {
    stop("'n' must be a positive whole number")
  }
  if (!vt_base_specified(object)) stop("the base process has unspecified coefficients; fit it to data first")
  as.numeric(vt_simulate(object, n, 1L))
})

#' @describeIn ktau Kendall's tau of the pair copulas of the base process of a
#'   [vtscop-class] model (the kpacf, at lags `1, ..., lagmax`, default 20).
#' @export
setMethod("ktau", "vtscop", function(object, lagmax, ...) {
  if (missing(lagmax)) lagmax <- 20
  ktau(object@base, lagmax)
})
