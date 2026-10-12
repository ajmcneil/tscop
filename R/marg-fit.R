# Maximum likelihood fitting of parametric margins.

#' Fitted parametric margins
#'
#' Class of objects returned by [tsfit()] for [marg-class] models.
#'
#' @slot model the fitted [marg-class] model.
#' @slot data the data to which the model was fitted (missing values removed).
#' @slot estimate named numeric vector of the estimates of the free parameters.
#' @slot se standard errors matching `estimate`; `NA` when not computed.
#' @slot vcov the matching covariance matrix; `NA` when not computed.
#' @slot se_method `"none"`, `"hessian"` or `"bootstrap"`.
#' @slot boot matrix of bootstrap estimates (no rows unless `se = "bootstrap"`).
#' @slot loglik maximised log-likelihood.
#' @slot nobs number of observations.
#' @slot npar number of estimated parameters.
#' @slot convergence convergence code of the optimisation (0 is success).
#' @slot message optimiser message.
#' @slot start the [marg-class] model with the starting values.
#' @slot template the original input if it was a `ts`, `zoo` or `xts` object,
#'   otherwise `NULL`.
#'
#' @exportClass margfit
setClass("margfit", slots = c(
  model = "marg", data = "numeric", estimate = "numeric", se = "numeric", vcov = "matrix",
  se_method = "character", boot = "matrix", loglik = "numeric", nobs = "integer",
  npar = "integer", convergence = "integer", message = "character", start = "marg",
  template = "ANY"
))

mg_negll <- function(z, fam, free, th_fixed, x) {
  if (any(!is.finite(z)) || any(abs(z) > 12)) return(1e10)
  th <- th_fixed
  th[free] <- mg_from_free(z, fam$types[free])
  if (th[["sigma"]] < 1e-300 || (!is.na(th["df"]) && th[["df"]] > 1e5)) return(1e10)
  ll <- sum(mg_ldens(fam, th, x))
  if (is.finite(ll)) -ll else 1e10
}

# Core ML fit of one family from a starting parameter vector; returns the full
# parameter vector, the free-scale optimum and the optimiser output.
mg_fit_core <- function(x, fam, th0, fix, control = list()) {
  free <- setdiff(fam$pars, fix)
  ctrl <- utils::modifyList(list(maxit = 500, reltol = 1e-10), control)
  z0 <- mg_to_free(th0[free], fam$types[free])
  if (length(free) == 0L) {
    return(list(th = th0, z = numeric(0), free = free, value = mg_negll(z0, fam, free, th0, x),
                convergence = 0L, message = ""))
  }
  if (mg_negll(z0, fam, free, th0, x) >= 1e10) stop("the starting values are infeasible for these data", call. = FALSE)
  c1 <- utils::modifyList(ctrl, list(ndeps = rep(1e-6, length(z0))))
  opt <- stats::optim(z0, mg_negll, fam = fam, free = free, th_fixed = th0, x = x, method = "BFGS", control = c1)
  if (opt$convergence != 0L || opt$value >= 1e10) {
    alt <- stats::optim(opt$par, mg_negll, fam = fam, free = free, th_fixed = th0, x = x, method = "Nelder-Mead",
                        control = list(maxit = 5000, reltol = 1e-10, warn.1d.NelderMead = FALSE))
    if (alt$value < opt$value) opt <- alt
  }
  th <- th0
  th[free] <- mg_from_free(opt$par, fam$types[free])
  list(th = th, z = opt$par, free = free, value = opt$value, convergence = opt$convergence,
       message = if (is.null(opt$message)) "" else opt$message)
}

# Starting values: the ladder fits the parent family first (with the shared
# fixed parameters), then maps its estimates to the child's.
mg_start <- function(x, fam, th_user, fix, control) {
  th0 <- fam$start(x)
  if (!is.na(fam$parent)) {
    pf <- mg_family(fam$parent)
    pth <- th_user[pf$pars]
    pth[is.na(pth)] <- NA
    pstart <- mg_start(x, pf, th_user, intersect(fix, pf$pars), control)
    pfit <- tryCatch(mg_fit_core(x, pf, pstart, intersect(fix, pf$pars), control), error = function(e) NULL)
    if (!is.null(pfit)) {
      th0 <- tryCatch(fam$from_parent(pfit$th), error = function(e) th0)
      th0 <- th0[fam$pars]
    }
  }
  # user values override for the parameters that were supplied
  ok <- !is.na(th_user)
  th0[ok] <- th_user[ok]
  th0
}

#' @describeIn tsfit Maximum likelihood fit of a [marg-class] model to the data
#'   `U` (here the raw observations, not copula data). The starting values come
#'   from `init` (a [marg-class] model or a named parameter vector), otherwise
#'   from the ladder: the nested parent family is fitted first (Gaussian for
#'   the Student t, Student t for the skew t, Laplace for the skew Laplace) and
#'   its estimates start the fit; parameters specified in `object` (and those in
#'   `fix`) are used as given. `se` is `"hessian"` or `"bootstrap"` (parametric
#'   bootstrap with refits from the fitted values), or `FALSE`.
#' @export
setMethod("tsfit", signature(U = "ANY", object = "marg"),
  function(U, object, se = FALSE, B = 200, init = NULL, control = list(), ...) {
    if (isTRUE(se)) se <- "bootstrap"
    if (!identical(se, FALSE) && !(is.character(se) && length(se) == 1L && se %in% c("hessian", "bootstrap"))) {
      stop("'se' must be FALSE, \"hessian\" or \"bootstrap\"", call. = FALSE)
    }
    if (identical(se, "bootstrap") && (!is.numeric(B) || length(B) != 1L || B < 10 || B != round(B))) {
      stop("'B' must be a whole number of at least 10", call. = FALSE)
    }
    if (!is.numeric(U)) stop("the data must be numeric", call. = FALSE)
    x <- as.numeric(U)
    x <- x[is.finite(x)]
    if (length(x) < 10L) stop("at least 10 finite observations are needed", call. = FALSE)
    fam <- mg_family(object@family)
    th_user <- object@pars
    if (!is.null(init)) {
      iv <- if (is(init, "marg")) init@pars else unlist(init)
      if (is.null(names(iv)) || any(!(names(iv) %in% fam$pars))) stop("'init' must be named with parameters of the family", call. = FALSE)
      th_user[names(iv)] <- iv
    }
    fix <- object@fix
    th_user[fix] <- object@pars[fix] # fixed parameters keep their values whatever 'init' says
    th0 <- if (!is.null(init) && !anyNA(th_user)) th_user else mg_start(x, fam, th_user, fix, control)
    if (anyNA(th0)) stop("could not determine starting values", call. = FALSE)
    fit <- mg_fit_core(x, fam, th0, fix, control)
    free <- fit$free
    model <- new("marg", family = object@family, pars = fit$th, fix = fix)
    est <- fit$th[free]
    npar <- length(free)
    V <- matrix(NA_real_, npar, npar, dimnames = list(free, free))
    boot <- matrix(NA_real_, 0L, npar, dimnames = list(NULL, free))
    se_method <- "none"
    negll <- function(z) mg_negll(z, fam, free, fit$th, x)
    to_nat <- function(z) mg_from_free(z, fam$types[free])
    if (identical(se, "hessian") && npar > 0L) {
      H <- num_hessian(negll, fit$z)
      Vz <- tryCatch(solve(H), error = function(e) NULL)
      if (is.null(Vz) || any(!is.finite(Vz)) || any(diag(Vz) <= 0) ||
          any(eigen(H, symmetric = TRUE, only.values = TRUE)$values <= 0)) {
        warning("the information matrix is not positive definite; standard errors not available")
      } else {
        J <- num_jacobian(to_nat, fit$z)
        Vn <- J %*% Vz %*% t(J)
        V[] <- (Vn + t(Vn)) / 2
        se_method <- "hessian"
      }
    } else if (identical(se, "bootstrap") && npar > 0L) {
      boot <- matrix(NA_real_, B, npar, dimnames = list(NULL, free))
      for (b in seq_len(B)) {
        xb <- mg_qf(fam, fit$th, stats::runif(length(x)))
        fb <- tryCatch(suppressWarnings(mg_fit_core(x = xb, fam = fam, th0 = fit$th, fix = fix, control = control)),
                       error = function(e) NULL)
        if (!is.null(fb) && fb$convergence == 0L) boot[b, ] <- fb$th[free]
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
    names(se_vec) <- free
    start_model <- new("marg", family = object@family, pars = th0, fix = fix)
    new("margfit",
      model = model, data = x, estimate = est, se = se_vec, vcov = V, se_method = se_method,
      boot = boot, loglik = -fit$value, nobs = length(x), npar = as.integer(npar),
      convergence = as.integer(fit$convergence), message = fit$message, start = start_model,
      template = if (inherits(U, c("ts", "zoo", "xts"))) U else NULL
    )
  }
)

# ---- methods for fitted margins ----------------------------------------------------

#' @rdname tscop-methods
#' @export
setMethod("coef", "margfit", function(object) object@estimate)

#' @rdname tscop-methods
#' @export
setMethod("logLik", "margfit", function(object) {
  structure(object@loglik, df = object@npar, nobs = object@nobs, class = "logLik")
})

#' @rdname tscop-methods
#' @export
setMethod("vcov", "margfit", function(object, ...) object@vcov)

#' @rdname tscop-methods
#' @export
setMethod("show", "margfit", function(object) {
  m <- object@model
  cat("object class: margfit\n")
  cat("family: ", m@family, " (", mg_family(m@family)$description, ")\n", sep = "")
  cat("estimates:\n")
  if (object@se_method == "none") {
    print(object@estimate)
  } else {
    print(cbind(estimate = object@estimate, se = object@se))
    cat("standard errors: ", object@se_method, "\n", sep = "")
  }
  if (length(m@fix) > 0L) cat("fixed: ", paste(paste0(m@fix, " = ", format(m@pars[m@fix])), collapse = ", "), "\n", sep = "")
  cat("loglik: ", format(object@loglik, digits = 7), "  npar: ", object@npar, "  nobs: ", object@nobs, "\n", sep = "")
  if (object@convergence != 0L) cat("warning: optimiser did not converge (code ", object@convergence, ")\n", sep = "")
  invisible(object)
})

#' @rdname tssim
#' @export
setMethod("tssim", "margfit", function(object, n, ...) {
  if (missing(n)) n <- length(object@data)
  tssim(object@model, n)
})

#' @rdname marg-functions
#' @export
setMethod("mdens", "margfit", function(object, x, log = FALSE, ...) mdens(object@model, x, log))
#' @rdname marg-functions
#' @export
setMethod("mdf", "margfit", function(object, x, lower.tail = TRUE, ...) mdf(object@model, x, lower.tail))
#' @rdname marg-functions
#' @export
setMethod("mqf", "margfit", function(object, p, lower.tail = TRUE, ...) mqf(object@model, p, lower.tail))

mg_fit_data <- function(object) if (is.null(object@template)) object@data else object@template

#' @rdname resid-method
#' @export
setMethod("resid", "margfit", function(object, data, scale = "uniform", nodes = 64, ...) {
  if (missing(data)) data <- mg_fit_data(object)
  mg_resid(object@model, data, scale)
})

#' @rdname resid-method
#' @export
setMethod("residuals", "margfit", function(object, data, scale = "uniform", nodes = 64, ...) {
  if (missing(data)) data <- mg_fit_data(object)
  mg_resid(object@model, data, scale)
})
