# Fitting and methods for arfactcop as a stand-alone Gaussian copula model.

#' Fitted AR-factor copula processes
#'
#' Class of objects returned by [tsfit()] for [arfactcop-class] models.
#'
#' @slot model the fitted [arfactcop-class] model.
#' @slot data the copula data to which the model was fitted.
#' @slot pseudo logical: were the data pseudo-observations?
#' @slot estimate named numeric vector of estimates.
#' @slot se standard errors matching `estimate`; `NA` when not computed.
#' @slot vcov the matching covariance matrix; `NA` when not computed.
#' @slot se_method `"none"`, `"hessian"` or `"bootstrap"`.
#' @slot boot matrix of bootstrap estimates (no rows unless
#'   `se = "bootstrap"`).
#' @slot loglik maximised copula log-likelihood.
#' @slot nobs number of non-missing observations.
#' @slot npar number of estimated parameters.
#' @slot convergence convergence code of the optimiser (0 is success).
#' @slot message optimiser message.
#' @slot start the [arfactcop-class] model used as starting values.
#' @slot template the original input if it was a `ts`, `zoo` or `xts` object,
#'   otherwise `NULL`.
#'
#' @exportClass arfactcopfit
setClass("arfactcopfit", slots = c(
  model = "arfactcop", data = "numeric", pseudo = "logical",
  estimate = "numeric", se = "numeric", vcov = "matrix", se_method = "character",
  boot = "matrix", loglik = "numeric", nobs = "integer", npar = "integer",
  convergence = "integer", message = "character", start = "arfactcop",
  template = "ANY"
))

# ---- simulation ---------------------------------------------------------

# Exact stationary simulation of the normal scores: m independent stationary
# unit-variance AR(1) processes plus white noise.
afc_simulate <- function(n, lambda, w, w0) {
  x <- sqrt(w0) * stats::rnorm(n)
  for (i in seq_along(lambda)) {
    e <- sqrt(1 - lambda[i]^2) * stats::rnorm(n)
    e[1L] <- stats::rnorm(1L)
    y <- as.numeric(stats::filter(e, lambda[i], method = "recursive"))
    x <- x + sqrt(w[i]) * y
  }
  x
}

afc_check_known <- function(model) {
  if (!gc_is_specified(model)) stop("the model has unspecified coefficients; fit it to data first")
  w0 <- if (model@noise) 1 - sum(model@w) else 0
  list(lambda = model@lambda, w = model@w, w0 = max(w0, 0))
}

#' @describeIn tssim Exact simulation from a known [arfactcop-class] model:
#'   independent stationary AR(1) components plus noise (`n` defaults to 1000).
#' @export
setMethod("tssim", "arfactcop", function(object, n, ...) {
  if (missing(n)) n <- 1000
  if (!is.numeric(n) || length(n) != 1L || n < 1 || n != round(n)) {
    stop("'n' must be a positive whole number")
  }
  p <- afc_check_known(object)
  stats::pnorm(afc_simulate(n, p$lambda, p$w, p$w0))
})

#' @describeIn ktau Kendall's tau of the Gaussian pair copulas of an
#'   [arfactcop-class] model at lags `1, ..., lagmax` (default 20).
#' @export
setMethod("ktau", "arfactcop", function(object, lagmax, ...) {
  if (missing(lagmax)) lagmax <- 20
  if (!is.numeric(lagmax) || length(lagmax) != 1L || lagmax < 1 || lagmax != round(lagmax)) {
    stop("'lagmax' must be a positive whole number")
  }
  gc_tau(object, lagmax)
})

# ---- residuals and prediction --------------------------------------------

afc_resid <- function(model, data, scale) {
  scale <- match.arg(scale, c("uniform", "normal"))
  u <- check_copula_data(data, min_n = 1L)
  p <- afc_check_known(model)
  f <- afc_filter(stats::qnorm(u), p$lambda, p$w, p$w0)
  r <- f$v / sqrt(f$F)
  restore_attrs(if (scale == "uniform") stats::pnorm(r) else r, data)
}

afc_predict <- function(model, data, h, type, x) {
  type <- gauss_predict_check(h, type, x)
  u <- check_copula_data(data, min_n = 1L)
  p <- afc_check_known(model)
  f <- afc_filter(stats::qnorm(u), p$lambda, p$w, p$w0)
  a <- f$a
  P <- f$P
  zv <- sqrt(p$w)
  ll2 <- outer(p$lambda, p$lambda)
  Q <- diag(1 - p$lambda^2, length(p$lambda))
  hmax <- max(h)
  mu <- numeric(hmax)
  v <- numeric(hmax)
  for (k in seq_len(hmax)) {
    mu[k] <- sum(zv * a)
    v[k] <- sum(zv * as.numeric(P %*% zv)) + p$w0
    a <- p$lambda * a
    P <- ll2 * P + Q
  }
  gauss_predict_result(mu, v, h, type, x)
}

#' @rdname resid-tscop
#' @export
setMethod("resid", "arfactcop", function(object, data, scale = "uniform", ...) {
  afc_resid(object, data, scale)
})

#' @rdname resid-tscop
#' @export
setMethod("residuals", "arfactcop", function(object, data, scale = "uniform", ...) {
  afc_resid(object, data, scale)
})

#' @rdname predict-tscop
#' @export
setMethod("predict", "arfactcop", function(object, data, h = 1, type = "moments",
                                           x = NULL, ...) {
  afc_predict(object, data, h, type, x)
})

# ---- fitting ---------------------------------------------------------------

#' @describeIn tsfit Maximum likelihood fit of an [arfactcop-class] model (exact
#'   Gaussian likelihood by a Kalman filter; `se` as for [armacop-class]
#'   models). Starting values come from `init` (an `arfactcop` model or a
#'   numeric vector in `coef` order), from `object` if its coefficients are
#'   specified, or from a least-squares fit to the sample autocorrelations of
#'   the normal scores.
#' @export
setMethod("tsfit", signature(U = "ANY", object = "arfactcop"),
  function(U, object, pseudo = TRUE, se = FALSE, B = 200, init = NULL,
           control = list(), ...) {
    if (isTRUE(se)) se <- "bootstrap"
    if (!identical(se, FALSE) && !(is.character(se) && length(se) == 1L &&
                                   se %in% c("hessian", "bootstrap"))) {
      stop("'se' must be FALSE, \"hessian\" or \"bootstrap\"")
    }
    if (identical(se, "bootstrap") &&
        (!is.numeric(B) || length(B) != 1L || B < 10 || B != round(B))) {
      stop("'B' must be a whole number of at least 10")
    }
    u <- check_copula_data(U)
    fit <- arfactcop_mle(u, object, init, control)
    model <- fit$model
    est <- coef(model)
    cn <- names(est)
    npar <- length(est)
    m <- object@m
    noise <- object@noise
    V <- matrix(NA_real_, npar, npar, dimnames = list(cn, cn))
    boot <- matrix(NA_real_, 0L, npar, dimnames = list(NULL, cn))
    se_method <- "none"
    if (identical(se, "hessian")) {
      ops <- gc_free_ops(object)
      H <- num_hessian(function(v) afc_negll(v, m, noise, fit$x), fit$z)
      Vz <- tryCatch(solve(H), error = function(e) NULL)
      if (is.null(Vz) || any(!is.finite(Vz)) || any(diag(Vz) <= 0) ||
          any(eigen(H, symmetric = TRUE, only.values = TRUE)$values <= 0)) {
        warning("the information matrix is not positive definite; standard errors not available")
      } else {
        J <- num_jacobian(function(v) ops$coef(v), fit$z)
        Vn <- J %*% Vz %*% t(J)
        V[] <- (Vn + t(Vn)) / 2
        se_method <- "hessian"
      }
    } else if (identical(se, "bootstrap")) {
      n <- length(u)
      missing <- is.na(u)
      boot <- matrix(NA_real_, B, npar, dimnames = list(NULL, cn))
      for (b in seq_len(B)) {
        ub <- tssim(model, n)
        ub[missing] <- NA
        if (pseudo) {
          ub <- rank(ub, na.last = "keep", ties.method = "random") / (sum(!missing) + 1)
        }
        fb <- tryCatch(suppressWarnings(arfactcop_mle(ub, object, init, control)$model),
                       error = function(e) NULL)
        if (!is.null(fb)) boot[b, ] <- coef(fb)
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
    names(se_vec) <- cn
    new("arfactcopfit",
      model = model, data = u, pseudo = pseudo, estimate = est, se = se_vec, vcov = V,
      se_method = se_method, boot = boot,
      loglik = -fit$opt$value - sum(stats::dnorm(fit$x, log = TRUE), na.rm = TRUE),
      nobs = as.integer(sum(!is.na(u))), npar = as.integer(npar),
      convergence = as.integer(fit$opt$convergence),
      message = if (is.null(fit$opt$message)) "" else fit$opt$message,
      start = fit$start,
      template = if (inherits(U, c("ts", "zoo", "xts"))) U else NULL
    )
  }
)

# ---- methods for fitted objects ---------------------------------------------

#' @rdname tscop-methods
#' @export
setMethod("coef", "arfactcopfit", function(object) object@estimate)

#' @rdname tscop-methods
#' @export
setMethod("logLik", "arfactcopfit", function(object) {
  structure(object@loglik, df = object@npar, nobs = object@nobs, class = "logLik")
})

#' @rdname tscop-methods
#' @export
setMethod("vcov", "arfactcopfit", function(object, ...) object@vcov)

#' @rdname tscop-methods
#' @export
setMethod("show", "arfactcopfit", function(object) {
  cat("object class: arfactcopfit\n")
  cat("name: ", object@model@name, "\n", sep = "")
  cat("estimates:\n")
  if (object@se_method == "none") {
    print(object@estimate)
  } else {
    print(cbind(estimate = object@estimate, se = object@se))
    cat("standard errors: ", object@se_method, "\n", sep = "")
  }
  if (object@model@noise) cat("noise share w0: ", format(1 - sum(object@model@w), digits = 4), "\n", sep = "")
  cat("loglik: ", format(object@loglik, digits = 7), "  npar: ", object@npar,
      "  nobs: ", object@nobs, "\n", sep = "")
  if (object@convergence != 0L) {
    cat("warning: optimiser did not converge (code ", object@convergence, ")\n", sep = "")
  }
  invisible(object)
})

#' @rdname tssim
#' @export
setMethod("tssim", "arfactcopfit", function(object, n, ...) {
  if (missing(n)) n <- length(object@data)
  tssim(object@model, n)
})

#' @rdname ktau
#' @export
setMethod("ktau", "arfactcopfit", function(object, lagmax, ...) {
  if (missing(lagmax)) lagmax <- 20
  ktau(object@model, lagmax)
})

afc_fit_data <- function(object) if (is.null(object@template)) object@data else object@template

#' @rdname resid-tscop
#' @export
setMethod("resid", "arfactcopfit", function(object, data, scale = "uniform", ...) {
  if (missing(data)) data <- afc_fit_data(object)
  afc_resid(object@model, data, scale)
})

#' @rdname resid-tscop
#' @export
setMethod("residuals", "arfactcopfit", function(object, data, scale = "uniform", ...) {
  if (missing(data)) data <- afc_fit_data(object)
  afc_resid(object@model, data, scale)
})

#' @rdname predict-tscop
#' @export
setMethod("predict", "arfactcopfit", function(object, data, h = 1, type = "moments",
                                              x = NULL, ...) {
  if (missing(data)) data <- object@data
  afc_predict(object@model, data, h, type, x)
})
