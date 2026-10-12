#' Parametric margins
#'
#' A parametric marginal distribution for a time series: a location-scale
#' family \eqn{x = \mu + \sigma z}, with `sigma` a pure scale parameter (not the
#' standard deviation) and shape parameters depending on the family. The
#' families are `"gauss"` (Gaussian), `"laplace"`, `"slaplace"` (skew
#' Laplace), `"st"` (Student t, `df` degrees of freedom), `"sst"` (skew
#' Student t), `"jsu"` (Johnson SU), `"nig"` (normal inverse Gaussian), `"hyp"`
#' (hyperbolic) and `"gh"` (generalised hyperbolic). The skew families use the Fernandez-Steel construction with the
#' skewing parameter `gamma` (`gamma = 1` is symmetric; `mu` is then the mode of
#' the distribution): the density of the standardised variable is
#' \eqn{2/(\gamma + 1/\gamma) f_0(\gamma z)} for \eqn{z < 0} and
#' \eqn{2/(\gamma + 1/\gamma) f_0(z/\gamma)} for \eqn{z \ge 0}, where \eqn{f_0}
#' is the symmetric base density. The Johnson SU distribution is defined by
#' \eqn{\gamma + \delta\,\mathrm{asinh}((x - \mu)/\sigma) \sim N(0, 1)}
#' with `delta` > 0 (smaller values give heavier tails; the Gaussian is the
#' limit as `delta` tends to infinity with `sigma` proportional to `delta`) and
#' `gamma` real (skewness; zero is symmetric about `mu`). The normal inverse
#' Gaussian distribution is \eqn{x = \mu + \sigma(\rho W + \sqrt{W} N)} with
#' \eqn{W} inverse Gaussian of mean one and shape `shape` and \eqn{N} standard
#' normal, and `skew` = \eqn{\rho} a scale-free skewness parameter
#' (`skew` = 0 is symmetric; the mean of the standardised variable is `skew`
#' and its variance is `1 + skew^2 / shape`; the Gaussian is the limit as `shape`
#' tends to infinity). This is the `alpha.bar` parameterisation of the `ghyp`
#' package with the skewness `gamma` of that package equal to `sigma * skew`,
#' so that `sigma` is a pure scale parameter. The generalised hyperbolic
#' distribution has the same representation with \eqn{W} generalised inverse
#' Gaussian of mean one, order `lambda` (any real number) and
#' shape parameter `shape` (its `alpha.bar`); `"nig"` is the case `lambda` = −1/2
#' and `"hyp"` the hyperbolic distribution, `lambda` = 1.
#'
#' Parameters named in `fix` are held fixed at their values in `pars` when the
#' model is fitted with [tsfit()]; for example
#' `marg("st", pars = c(mu = 0), fix = "mu")` is a zero-mean Student t.
#'
#' @slot family name of the family.
#' @slot pars named numeric vector of all the parameters of the family
#'   (`mu`, `sigma`, then the shape parameters); `NA` for parameters not yet
#'   specified (a model structure to be fitted).
#' @slot fix names of the parameters held fixed in fitting.
#'
#' @exportClass marg
setClass("marg", slots = c(family = "character", pars = "numeric", fix = "character"))

mg_bounds_ok <- function(fam, pars) {
  th <- pars[!is.na(pars)]
  if ("sigma" %in% names(th) && th[["sigma"]] <= 0) return("'sigma' must be positive")
  if ("df" %in% names(th) && th[["df"]] <= 1) return("'df' must exceed 1")
  if ("gamma" %in% names(th) && fam$types[["gamma"]] == "pos" && th[["gamma"]] <= 0) return("'gamma' must be positive")
  if ("delta" %in% names(th) && th[["delta"]] <= 0) return("'delta' must be positive")
  if ("shape" %in% names(th) && th[["shape"]] <= 0) return("'shape' must be positive")
  if (any(!is.finite(th))) return("parameters must be finite")
  NULL
}

setValidity("marg", function(object) {
  fam <- .marg_families[[object@family]]
  if (is.null(fam)) return(paste0("unknown family '", object@family, "'"))
  if (!identical(names(object@pars), fam$pars)) {
    return(paste0("'pars' must be named ", paste(fam$pars, collapse = ", ")))
  }
  if (!all(object@fix %in% fam$pars)) return("'fix' must name parameters of the family")
  if (any(is.na(object@pars[object@fix]))) return("fixed parameters need values in 'pars'")
  msg <- mg_bounds_ok(fam, object@pars)
  if (!is.null(msg)) return(msg)
  TRUE
})

#' Constructor for parametric margins
#'
#' @param family one of `"gauss"`, `"laplace"`, `"slaplace"`, `"st"`, `"sst"`,
#'   `"jsu"`, `"nig"`, `"hyp"`, `"gh"`.
#' @param pars named numeric vector (or list) of parameter values; parameters
#'   that are omitted are unspecified (`NA`), giving a model structure to be
#'   fitted. See [marg-class] for the parameters of each family.
#' @param fix character vector of parameters held fixed (at their values in
#'   `pars`) when the model is fitted.
#'
#' @return An object of class [marg-class].
#' @export
#' @examples
#' marg("sst", c(mu = 0, sigma = 1, df = 5, gamma = 0.9))
#' marg("st")                                       # structure to be fitted
#' marg("st", c(mu = 0), fix = "mu")                # zero-mean Student t
marg <- function(family, pars = NULL, fix = character(0)) {
  fam <- mg_family(family)
  th <- stats::setNames(rep(NA_real_, length(fam$pars)), fam$pars)
  if (!is.null(pars)) {
    pars <- unlist(pars)
    if (is.null(names(pars)) || any(!(names(pars) %in% fam$pars))) {
      stop("'pars' must be named with parameters of the family: ", paste(fam$pars, collapse = ", "), call. = FALSE)
    }
    th[names(pars)] <- as.numeric(pars)
  }
  new("marg", family = family, pars = th, fix = as.character(fix))
}

mg_specified <- function(object) !anyNA(object@pars)

mg_free_names <- function(object) setdiff(names(object@pars), object@fix)

mg_check_specified <- function(object) {
  if (!mg_specified(object)) stop("the margin has unspecified parameters; fit it to data first", call. = FALSE)
}

#' Density, distribution function and quantile function of a margin
#'
#' Generics evaluating the density (`mdens`), distribution function (`mdf`) and
#' quantile function (`mqf`) of a [marg-class] model with known parameters, or
#' of a fitted margin.
#'
#' @param object a [marg-class] model with specified parameters or a
#'   [margfit-class] object.
#' @param x values at which to evaluate the density or distribution function.
#' @param p probabilities (or upper-tail probabilities if `lower.tail = FALSE`).
#' @param log logical; return the log-density?
#' @param lower.tail logical; if `FALSE` the upper tail is used (accurate for
#'   extreme values).
#' @param ... not used.
#'
#' @return A numeric vector.
#' @name marg-functions
NULL

#' @rdname marg-functions
#' @export
setGeneric("mdens", function(object, x, log = FALSE, ...) standardGeneric("mdens"), signature = "object")
#' @rdname marg-functions
#' @export
setGeneric("mdf", function(object, x, lower.tail = TRUE, ...) standardGeneric("mdf"), signature = "object")
#' @rdname marg-functions
#' @export
setGeneric("mqf", function(object, p, lower.tail = TRUE, ...) standardGeneric("mqf"), signature = "object")

#' @rdname marg-functions
#' @export
setMethod("mdens", "marg", function(object, x, log = FALSE, ...) {
  mg_check_specified(object)
  out <- mg_ldens(mg_family(object@family), object@pars, as.numeric(x))
  if (log) out else exp(out)
})

#' @rdname marg-functions
#' @export
setMethod("mdf", "marg", function(object, x, lower.tail = TRUE, ...) {
  mg_check_specified(object)
  mg_cdf(mg_family(object@family), object@pars, as.numeric(x), lower.tail)
})

#' @rdname marg-functions
#' @export
setMethod("mqf", "marg", function(object, p, lower.tail = TRUE, ...) {
  mg_check_specified(object)
  p <- as.numeric(p)
  if (any(p < 0 | p > 1, na.rm = TRUE)) stop("'p' must be in [0, 1]", call. = FALSE)
  mg_qf(mg_family(object@family), object@pars, p, lower.tail)
})

#' @rdname tscop-methods
#' @export
setMethod("coef", "marg", function(object) object@pars[mg_free_names(object)])

#' @rdname tscop-methods
#' @export
setMethod("show", "marg", function(object) {
  fam <- mg_family(object@family)
  cat("object class: marg\n")
  cat("family: ", object@family, " (", fam$description, ")\n", sep = "")
  cat("parameters:\n")
  print(object@pars)
  if (length(object@fix) > 0L) cat("fixed: ", paste(object@fix, collapse = ", "), "\n", sep = "")
  invisible(object)
})

#' @describeIn tssim Random draws from a [marg-class] model with known
#'   parameters (`n` defaults to 1000).
#' @export
setMethod("tssim", "marg", function(object, n, ...) {
  if (missing(n)) n <- 1000
  if (!is.numeric(n) || length(n) != 1L || n < 1 || n != round(n)) stop("'n' must be a positive whole number")
  mg_check_specified(object)
  mg_rand(mg_family(object@family), object@pars, n)
})

mg_resid <- function(model, data, scale) {
  scale <- match.arg(scale, c("uniform", "normal"))
  mg_check_specified(model)
  if (!is.numeric(data)) stop("'data' must be numeric", call. = FALSE)
  u <- mg_cdf(mg_family(model@family), model@pars, as.numeric(data))
  u <- pmin(pmax(u, 1e-12), 1 - 1e-12)
  restore_attrs(if (scale == "uniform") u else stats::qnorm(u), data)
}

#' @rdname resid-method
#' @export
setMethod("resid", "marg", function(object, data, scale = "uniform", nodes = 64, ...) {
  mg_resid(object, data, scale)
})

#' @rdname resid-method
#' @export
setMethod("residuals", "marg", function(object, data, scale = "uniform", nodes = 64, ...) {
  mg_resid(object, data, scale)
})
