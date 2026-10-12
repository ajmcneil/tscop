#' Full time series models: a margin and a copula process
#'
#' A model for a time series \eqn{x_1, \ldots, x_n} with marginal distribution
#' \eqn{F} (a [marg-class] parametric margin or an [emarg-class] empirical margin)
#' and a copula process for the serial dependence of the probability integral
#' transforms \eqn{u_t = F(x_t)} (an [armacop-class], [arfactcop-class],
#' [sdvinecop-class] or [vtscop-class] model). The density of the data is
#' \eqn{c(u_1, \ldots, u_n) \prod_t f(x_t)}.
#'
#' Simulation applies the quantile function of the margin to a simulated copula
#' series. Estimation by [tsfit()] is in two stages: the parametric margin is
#' fitted by maximum likelihood and the copula process to its probability
#' integral transforms; with an empirical margin the copula process is fitted to
#' the pseudo-observations `rank / (n + 1)` (the semiparametric approach of
#' Bladt, Dias, Han and McNeil, 2026).
#'
#' @slot margin the margin: a [marg-class] or [emarg-class] object.
#' @slot copula the copula process: an [armacop-class], [arfactcop-class],
#'   [sdvinecop-class] or [vtscop-class] object.
#'
#' @references Bladt, M., Dias, A., Han, J. and McNeil, A. J. (2026).
#'   Semiparametric and parametric distributional forecasting of univariate time
#'   series using non-Gaussian ARMA models based on D-vines. *Canadian Journal
#'   of Statistics*, 54(1), e70033. \doi{10.1002/cjs.70033}
#'
#' @exportClass tsmod
setClass("tsmod", slots = c(margin = "ANY", copula = "ANY"))

tm_copula_ok <- function(cop) is(cop, "gausscop") || is(cop, "sdvinecop") || is(cop, "vtscop")

setValidity("tsmod", function(object) {
  if (!(is(object@margin, "marg") || is(object@margin, "emarg"))) {
    return("'margin' must be a marg or emarg object")
  }
  if (!tm_copula_ok(object@copula)) {
    return("'copula' must be an armacop, arfactcop, sdvinecop or vtscop object")
  }
  TRUE
})

#' Constructor for full time series models
#'
#' @param margin a [marg-class] or [emarg-class] object, or the name of a
#'   parametric family (`"sst"` is the same as `marg("sst")`, a margin to be
#'   fitted).
#' @param copula an [armacop-class], [arfactcop-class], [sdvinecop-class] or
#'   [vtscop-class] model, with known coefficients or a structure to be fitted.
#'
#' @return An object of class [tsmod-class].
#' @export
#' @examples
#' m <- tsmod(marg("sst", c(mu = 0, sigma = 1, df = 6, gamma = 1.2)),
#'            arfactcop(lambda = 0.9, w = 0.6))
#' m
#' set.seed(1)
#' x <- tssim(m, 500)
#' tsmod("sst", arfactcop(order = 1))      # a structure to be fitted
tsmod <- function(margin, copula) {
  if (is.character(margin) && length(margin) == 1L) margin <- marg(margin)
  new("tsmod", margin = margin, copula = copula)
}

#' @rdname tscop-methods
#' @export
setMethod("show", "tsmod", function(object) {
  cat("object class: tsmod\n")
  cat("margin:\n")
  cat(paste0("  ", utils::capture.output(show(object@margin))), sep = "\n")
  cat("copula process:\n")
  cat(paste0("  ", utils::capture.output(show(object@copula))), sep = "\n")
  invisible(object)
})

#' @rdname tscop-methods
#' @export
setMethod("coef", "tsmod", function(object) {
  cm <- if (is(object@margin, "marg")) coef(object@margin) else numeric(0)
  cc <- coef(object@copula)
  c(stats::setNames(cm, if (length(cm)) paste0("margin.", names(cm))),
    stats::setNames(cc, if (length(cc)) paste0("copula.", names(cc))))
})

#' @describeIn tssim Simulation from a [tsmod-class] model with known margin and
#'   copula: the quantile function of the margin applied to a simulated copula
#'   series.
#' @export
setMethod("tssim", "tsmod", function(object, n, ...) {
  if (missing(n)) n <- 1000
  u <- tssim(object@copula, n)
  x <- mqf(object@margin, u)
  x
})

#' @describeIn ktau Kendall's tau of the pair copulas of the copula process of a
#'   [tsmod-class] model.
#' @export
setMethod("ktau", "tsmod", function(object, lagmax, ...) {
  if (missing(lagmax)) lagmax <- 20
  ktau(object@copula, lagmax)
})

# Probability integral transforms of data under a margin (clamped to the open
# unit interval), and the pseudo-observations of an empirical margin.
tm_pits <- function(margin, x) {
  u <- mdf(margin, as.numeric(x))
  pmin(pmax(u, 1e-12), 1 - 1e-12)
}

tm_check_data <- function(x) {
  if (!is.numeric(x)) stop("the data must be numeric", call. = FALSE)
  if (anyNA(x) || any(!is.finite(x))) stop("the data must be complete and finite", call. = FALSE)
  if (length(x) < 20L) stop("at least 20 observations are needed", call. = FALSE)
  as.numeric(x)
}

#' @rdname resid-method
#' @param type for a [tsmod-class] model: `"copula"` (default), the residuals of
#'   the copula process applied to the probability integral transforms of the
#'   data, which are iid uniform if the model is correct, or `"margin"`, the
#'   probability integral transforms themselves.
#' @export
setMethod("resid", "tsmod", function(object, data, scale = "uniform", nodes = 64, type = "copula", ...) {
  tm_resid(object, data, scale, nodes, type)
})

#' @rdname resid-method
#' @export
setMethod("residuals", "tsmod", function(object, data, scale = "uniform", nodes = 64, type = "copula", ...) {
  tm_resid(object, data, scale, nodes, type)
})

tm_resid <- function(model, data, scale, nodes, type) {
  type <- match.arg(type, c("copula", "margin"))
  scale <- match.arg(scale, c("uniform", "normal"))
  if (!is.numeric(data)) stop("'data' must be numeric", call. = FALSE)
  u <- tm_pits(model@margin, data)
  out <- if (type == "margin") {
    if (scale == "uniform") u else stats::qnorm(u)
  } else {
    resid(model@copula, u, scale = scale, nodes = nodes)
  }
  restore_attrs(as.numeric(out), data)
}

#' Fitted full time series models
#'
#' Class of objects returned by [tsfit()] for [tsmod-class] models.
#'
#' @slot model the fitted [tsmod-class] model.
#' @slot data the data (original scale).
#' @slot margfit the fitted margin ([margfit-class]; the [emarg-class] object itself
#'   for an empirical margin).
#' @slot copfit the fitted copula process (its fitted class).
#' @slot estimate named numeric vector of the estimates, with the prefixes
#'   `margin.` and `copula.`.
#' @slot se,vcov,se_method standard errors, covariance matrix and their method
#'   (`NA` and `"none"` when not computed).
#' @slot boot matrix of bootstrap estimates.
#' @slot loglik the log-likelihood: for a parametric margin the sum of the margin
#'   and copula log-likelihoods, for an empirical margin the copula (pseudo)
#'   log-likelihood only.
#' @slot nobs number of observations.
#' @slot npar number of estimated parameters (both stages).
#' @slot method estimation method (`"twostage"`).
#' @slot template the original input if it was a `ts`, `zoo` or `xts` object,
#'   otherwise `NULL`.
#'
#' @exportClass tsmodfit
setClass("tsmodfit", slots = c(
  model = "tsmod", data = "numeric", margfit = "ANY", copfit = "ANY", estimate = "numeric",
  se = "numeric", vcov = "matrix", se_method = "character", boot = "matrix", loglik = "numeric",
  nobs = "integer", npar = "integer", method = "character", template = "ANY"
))

#' @describeIn tsfit Estimation of a [tsmod-class] model from the data `U` on its
#'   original scale (a numeric vector or time series without missing values).
#'   With a parametric margin the estimation is in two stages: the margin is
#'   fitted by maximum likelihood to the data and the copula process to the
#'   probability integral transforms (`tsfit(u, copula, pseudo = FALSE)`); with an
#'   empirical margin the copula process is fitted to the pseudo-observations
#'   (the semiparametric approach). Further arguments in `...` (for example
#'   `control` or `init`) are passed to the copula's fit; `margin_control` is the
#'   control list of the margin's optimiser. Standard errors are not yet
#'   available (`se = FALSE`).
#' @param method estimation method for a [tsmod-class] model; only `"twostage"`
#'   so far.
#' @param margin_control control list passed to the optimiser of the margin.
#' @export
setMethod("tsfit", signature(U = "ANY", object = "tsmod"),
  function(U, object, method = "twostage", se = FALSE, margin_control = list(), ...) {
    method <- match.arg(method, c("twostage"))
    if (!identical(se, FALSE)) stop("standard errors for tsmod models are not yet available (use se = FALSE)", call. = FALSE)
    x <- tm_check_data(U)
    template <- if (inherits(U, c("ts", "zoo", "xts"))) U else NULL
    mg <- object@margin
    if (is(mg, "marg")) {
      mfit <- tsfit(x, mg, control = margin_control)
      mmodel <- mfit@model
      u <- tm_pits(mmodel, x)
    } else {
      mfit <- emarg(x, type = mg@type, bw = if (mg@bwsel == "numeric") mg@bw else mg@bwsel,
                    qtype = mg@qtype, proper = mg@proper)
      mmodel <- mfit
      u <- tm_pits(mfit, x)
    }
    cfit <- tsfit(u, object@copula, pseudo = FALSE, ...)
    cmodel <- cfit@model
    model <- new("tsmod", margin = mmodel, copula = cmodel)
    cm <- if (is(mfit, "margfit")) coef(mfit) else numeric(0)
    cc <- coef(cfit)
    est <- c(stats::setNames(cm, if (length(cm)) paste0("margin.", names(cm))),
             stats::setNames(cc, if (length(cc)) paste0("copula.", names(cc))))
    npar <- length(est)
    ll <- as.numeric(logLik(cfit)) + if (is(mfit, "margfit")) as.numeric(logLik(mfit)) else 0
    new("tsmodfit",
      model = model, data = x, margfit = mfit, copfit = cfit, estimate = est,
      se = rep(NA_real_, npar), vcov = matrix(NA_real_, npar, npar, dimnames = list(names(est), names(est))),
      se_method = "none", boot = matrix(NA_real_, 0L, npar, dimnames = list(NULL, names(est))),
      loglik = ll, nobs = length(x), npar = as.integer(npar), method = method, template = template
    )
  }
)

#' @rdname tscop-methods
#' @export
setMethod("coef", "tsmodfit", function(object) object@estimate)

#' @rdname tscop-methods
#' @export
setMethod("logLik", "tsmodfit", function(object) {
  structure(object@loglik, df = object@npar, nobs = object@nobs, class = "logLik")
})

#' @rdname tscop-methods
#' @export
setMethod("vcov", "tsmodfit", function(object, ...) object@vcov)

#' @rdname tscop-methods
#' @export
setMethod("show", "tsmodfit", function(object) {
  cat("object class: tsmodfit   (method: ", object@method, ")\n", sep = "")
  cat("margin:\n")
  cat(paste0("  ", utils::capture.output(show(object@margfit))), sep = "\n")
  cat("copula process:\n")
  cat(paste0("  ", utils::capture.output(show(object@copfit))), sep = "\n")
  cat("loglik: ", format(object@loglik, digits = 7), "  npar: ", object@npar, "  nobs: ", object@nobs,
      if (!is(object@margfit, "margfit")) "  (copula pseudo-likelihood)", "\n", sep = "")
  invisible(object)
})

tm_fit_data <- function(object) if (is.null(object@template)) object@data else object@template

#' @rdname tssim
#' @export
setMethod("tssim", "tsmodfit", function(object, n, ...) {
  if (missing(n)) n <- length(object@data)
  tssim(object@model, n)
})

#' @rdname ktau
#' @export
setMethod("ktau", "tsmodfit", function(object, lagmax, ...) {
  if (missing(lagmax)) lagmax <- 20
  ktau(object@model, lagmax)
})

#' @rdname resid-method
#' @export
setMethod("resid", "tsmodfit", function(object, data, scale = "uniform", nodes = 64, type = "copula", ...) {
  if (missing(data)) data <- tm_fit_data(object)
  tm_resid(object@model, data, scale, nodes, type)
})

#' @rdname resid-method
#' @export
setMethod("residuals", "tsmodfit", function(object, data, scale = "uniform", nodes = 64, type = "copula", ...) {
  if (missing(data)) data <- tm_fit_data(object)
  tm_resid(object@model, data, scale, nodes, type)
})
