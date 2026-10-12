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

tm_auto <- function(margin) is.character(margin) && identical(margin, "auto")

setValidity("tsmod", function(object) {
  if (!(is(object@margin, "marg") || is(object@margin, "emarg") || tm_auto(object@margin))) {
    return("'margin' must be a marg or emarg object (or \"auto\")")
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
#'   fitted), or `"auto"`: the margin is then chosen from the data by
#'   [automarg()] when the model is fitted with [tsfit()] (the selection is
#'   recorded in the fitted object).
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
  if (is.character(margin) && length(margin) == 1L && !tm_auto(margin)) margin <- marg(margin)
  new("tsmod", margin = margin, copula = copula)
}

#' @rdname tscop-methods
#' @export
setMethod("show", "tsmod", function(object) {
  cat("object class: tsmod\n")
  cat("margin:\n")
  if (tm_auto(object@margin)) {
    cat("  chosen from the data by automarg when the model is fitted\n")
  } else {
    cat(paste0("  ", utils::capture.output(show(object@margin))), sep = "\n")
  }
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
  if (tm_auto(object@margin)) stop("the margin is chosen automatically when the model is fitted; fit it to data first", call. = FALSE)
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
#' @slot selection the [automarg-class] object if the margin was chosen
#'   automatically, otherwise `NULL`.
#' @slot template the original input if it was a `ts`, `zoo` or `xts` object,
#'   otherwise `NULL`.
#'
#' @exportClass tsmodfit
setClass("tsmodfit", slots = c(
  model = "tsmod", data = "numeric", margfit = "ANY", copfit = "ANY", estimate = "numeric",
  se = "numeric", vcov = "matrix", se_method = "character", boot = "matrix", loglik = "numeric",
  nobs = "integer", npar = "integer", method = "character", template = "ANY", selection = "ANY"
))

#' @describeIn tsfit Estimation of a [tsmod-class] model from the data `U` on its
#'   original scale (a numeric vector or time series without missing values).
#'   With a parametric margin the estimation is in two stages: the margin is
#'   fitted by maximum likelihood to the data and the copula process to the
#'   probability integral transforms (`tsfit(u, copula, pseudo = FALSE)`); with an
#'   empirical margin the copula process is fitted to the pseudo-observations
#'   (the semiparametric approach); with `margin = "auto"` the parametric margin
#'   is first chosen by [automarg()] (options in `automarg_args`, for example
#'   `list(families = c("st", "sst", "nig"), criterion = "BIC")`) and the
#'   selection is stored in the fitted object. Further arguments in `...` (for
#'   example `control` or `init`) are passed to the copula's fit;
#'   `margin_control` is the control list of the margin's optimiser. Standard
#'   errors come from the parametric bootstrap (`se = "bootstrap"`): `B` series
#'   are simulated from the fitted model and both stages are refitted from the
#'   fitted values, so that the estimation of the margin enters the standard
#'   errors of the copula parameters; with an automatically chosen margin the
#'   selected family is kept in the refits. With `cores > 1` the refits run in
#'   parallel (`parallel::mclapply`, not on Windows); the simulated series are
#'   generated first, so the results do not depend on `cores`.
#' @param method estimation method for a [tsmod-class] model; only `"twostage"`
#'   so far.
#' @param margin_control control list passed to the optimiser of the margin.
#' @param automarg_args list of arguments of [automarg()] (`families`, `fix`,
#'   `criterion`) used when the margin of a [tsmod-class] model is `"auto"`.
#' @param cores number of processes for the bootstrap refits of a
#'   [tsmod-class] model (default 1).
#' @export
setMethod("tsfit", signature(U = "ANY", object = "tsmod"),
  function(U, object, method = "twostage", se = FALSE, B = 200, cores = 1L,
           margin_control = list(), automarg_args = list(), ...) {
    method <- match.arg(method, c("twostage"))
    if (isTRUE(se)) se <- "bootstrap"
    if (!identical(se, FALSE) && !identical(se, "bootstrap")) {
      stop("'se' must be FALSE or \"bootstrap\" (Hessian standard errors are not available for two-stage fits)", call. = FALSE)
    }
    if (identical(se, "bootstrap") && (!is.numeric(B) || length(B) != 1L || B < 10 || B != round(B))) {
      stop("'B' must be a whole number of at least 10", call. = FALSE)
    }
    x <- tm_check_data(U)
    template <- if (inherits(U, c("ts", "zoo", "xts"))) U else NULL
    mg <- object@margin
    selection <- NULL
    if (tm_auto(mg)) {
      selection <- do.call(automarg, c(list(x = x, control = margin_control), automarg_args))
      mfit <- selection@best
      mmodel <- mfit@model
      u <- tm_pits(mmodel, x)
    } else if (is(mg, "marg")) {
      mfit <- tsfit(x, mg, control = margin_control)
      mmodel <- mfit@model
      u <- tm_pits(mmodel, x)
    } else {
      mfit <- tm_emarg_like(mg, x)
      mmodel <- mfit
      u <- tm_pits(mfit, x)
    }
    cfit <- tsfit(u, object@copula, pseudo = FALSE, ...)
    model <- new("tsmod", margin = mmodel, copula = cfit@model)
    est <- tm_coef(mfit, cfit)
    npar <- length(est)
    ll <- as.numeric(logLik(cfit)) + if (is(mfit, "margfit")) as.numeric(logLik(mfit)) else 0
    V <- matrix(NA_real_, npar, npar, dimnames = list(names(est), names(est)))
    boot <- matrix(NA_real_, 0L, npar, dimnames = list(NULL, names(est)))
    se_method <- "none"
    if (identical(se, "bootstrap") && npar > 0L) {
      n <- length(x)
      Xb <- lapply(seq_len(B), function(b) tssim(model, n))     # simulated in the main process
      refit <- function(xb) {
        tryCatch(suppressWarnings({
          mfb <- if (is(mmodel, "marg")) tsfit(xb, mmodel, control = margin_control) else tm_emarg_like(mmodel, xb)
          mm <- if (is(mfb, "margfit")) mfb@model else mfb
          cfb <- tsfit(tm_pits(mm, xb), model@copula, pseudo = FALSE, ...)
          tm_coef(mfb, cfb)
        }), error = function(e) NULL)
      }
      res <- if (cores > 1L && .Platform$OS.type != "windows") {
        parallel::mclapply(Xb, refit, mc.cores = cores)
      } else {
        lapply(Xb, refit)
      }
      boot <- matrix(NA_real_, B, npar, dimnames = list(NULL, names(est)))
      for (b in seq_len(B)) {
        if (is.numeric(res[[b]]) && length(res[[b]]) == npar) boot[b, ] <- res[[b]]
      }
      good <- stats::complete.cases(boot)
      if (sum(good) >= 2L) {
        V[] <- stats::cov(boot[good, , drop = FALSE])
        se_method <- "bootstrap"
        if (sum(good) < B) warning(sprintf("%d of %d bootstrap refits failed", B - sum(good), B))
      } else {
        warning("too few successful bootstrap refits; standard errors not available")
      }
    }
    se_vec <- if (se_method == "none") rep(NA_real_, npar) else sqrt(diag(V))
    names(se_vec) <- names(est)
    new("tsmodfit",
      model = model, data = x, margfit = mfit, copfit = cfit, estimate = est,
      se = se_vec, vcov = V, se_method = se_method, boot = boot,
      loglik = ll, nobs = length(x), npar = as.integer(npar), method = method, template = template,
      selection = selection
    )
  }
)

# an empirical margin with the settings of `template`, re-estimated on x
tm_emarg_like <- function(template, x) {
  emarg(x, type = template@type, bw = if (template@bwsel == "numeric") template@bw else template@bwsel,
        qtype = template@qtype, proper = template@proper)
}

# combined named coefficient vector of the stage fits
tm_coef <- function(mfit, cfit) {
  cm <- if (is(mfit, "margfit")) coef(mfit) else numeric(0)
  cc <- coef(cfit)
  c(stats::setNames(cm, if (length(cm)) paste0("margin.", names(cm))),
    stats::setNames(cc, if (length(cc)) paste0("copula.", names(cc))))
}

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
  if (!is.null(object@selection)) {
    cat("margin chosen by automarg (", object@selection@criterion, "): ", object@model@margin@family,
        "; ", nrow(object@selection@table), " families compared\n", sep = "")
  }
  cat("margin:\n")
  cat(paste0("  ", utils::capture.output(show(object@margfit))), sep = "\n")
  cat("copula process:\n")
  cat(paste0("  ", utils::capture.output(show(object@copfit))), sep = "\n")
  if (object@se_method != "none") {
    cat("combined estimates (standard errors: ", object@se_method, "):\n", sep = "")
    print(cbind(estimate = object@estimate, se = object@se))
  }
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
