#' Automatic selection of a parametric margin
#'
#' Fits several parametric margin families to the data and ranks them by AIC or
#' BIC. The families are fitted in order of nesting, and each fit starts from the
#' estimates of its nested parent family when that is among the families
#' (Gaussian -> Student t -> skew Student t; Gaussian -> Johnson SU; Laplace -> skew Laplace), which makes
#' the larger families easier to estimate. Families that fail to fit are
#' reported, not fatal.
#'
#' @param x numeric data (the raw observations; non-finite values are
#'   dropped).
#' @param families character vector of families to try; by default all the
#'   available families (`"gauss"`, `"laplace"`, `"slaplace"`, `"st"`, `"sst"`, `"jsu"`).
#' @param fix a named numeric vector of parameter values held fixed in every
#'   family that has the parameter, for example `c(mu = 0)` for zero-mean
#'   margins.
#' @param criterion `"AIC"` (default) or `"BIC"`.
#' @param se standard errors for the selected model only: `FALSE` (default),
#'   `"hessian"` or `"bootstrap"`.
#' @param B number of bootstrap replicates for `se = "bootstrap"`.
#' @param control list of control parameters passed to the optimiser.
#'
#' @return An object of class [automarg-class]: the selected [margfit-class] model
#'   (`@best`), a ranking table of the families (`@table`: log-likelihood, number
#'   of parameters, AIC, BIC, the difference from the best criterion value and the
#'   convergence code; failures show `NA`), all the fits (`@fits`) and the
#'   criterion used.
#' @export
#' @examples
#' set.seed(1)
#' x <- tssim(marg("st", c(mu = 0, sigma = 1, df = 4)), 1000)
#' automarg(x)
automarg <- function(x, families = names(.marg_families), fix = NULL, criterion = c("AIC", "BIC"),
                     se = FALSE, B = 200, control = list()) {
  criterion <- match.arg(criterion)
  if (!is.numeric(x)) stop("the data must be numeric", call. = FALSE)
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) < 10L) stop("at least 10 finite observations are needed", call. = FALSE)
  families <- unique(as.character(families))
  if (length(families) < 1L) stop("'families' must name at least one family", call. = FALSE)
  for (f in families) mg_family(f) # errors for unknown names
  if (!is.null(fix)) {
    if (!is.numeric(fix) || is.null(names(fix)) || any(!(names(fix) %in% c("mu", "sigma", "df", "gamma", "delta")))) {
      stop("'fix' must be a named numeric vector of parameter values", call. = FALSE)
    }
  }
  # parents first
  depth <- function(f) {
    d <- 0L
    while (!is.na(par <- .marg_families[[f]]$parent)) {
      d <- d + 1L
      f <- par
    }
    d
  }
  families <- families[order(vapply(families, depth, 0L))]
  fits <- vector("list", length(families))
  names(fits) <- families
  notes <- stats::setNames(rep("", length(families)), families)
  for (f in families) {
    fam <- mg_family(f)
    fx <- if (is.null(fix)) fix else fix[names(fix) %in% fam$pars]
    obj <- marg(f, pars = if (length(fx)) fx else NULL, fix = names(fx))
    init <- NULL
    pf <- fam$parent
    if (!is.na(pf) && !is.null(fits[[pf]])) {
      init <- tryCatch(fam$from_parent(fits[[pf]]@model@pars)[fam$pars], error = function(e) NULL)
    }
    res <- tryCatch(suppressWarnings(tsfit(x, obj, init = init, control = control)),
                    error = function(e) e)
    if (inherits(res, "error")) {
      notes[f] <- conditionMessage(res)
    } else {
      fits[[f]] <- res
    }
  }
  tab <- data.frame(
    family = families,
    npar = vapply(families, function(f) if (is.null(fits[[f]])) NA_integer_ else fits[[f]]@npar, 0L),
    loglik = vapply(families, function(f) if (is.null(fits[[f]])) NA_real_ else fits[[f]]@loglik, 0),
    convergence = vapply(families, function(f) if (is.null(fits[[f]])) NA_integer_ else fits[[f]]@convergence, 0L),
    stringsAsFactors = FALSE, row.names = NULL
  )
  tab$AIC <- -2 * tab$loglik + 2 * tab$npar
  tab$BIC <- -2 * tab$loglik + log(length(x)) * tab$npar
  crit <- tab[[criterion]]
  if (all(is.na(crit))) stop("no family could be fitted: ", paste(unique(notes[notes != ""]), collapse = "; "), call. = FALSE)
  tab$delta <- crit - min(crit, na.rm = TRUE)
  tab$note <- unname(notes[families])
  tab <- tab[order(tab$delta, na.last = TRUE), ]
  rownames(tab) <- NULL
  bestname <- tab$family[1L]
  best <- fits[[bestname]]
  if (!identical(se, FALSE)) {
    best <- tsfit(x, marg(bestname, best@model@pars[best@model@fix], fix = best@model@fix),
                  init = best@model, se = se, B = B, control = control)
  }
  new("automarg", best = best, table = tab, fits = fits[!vapply(fits, is.null, TRUE)],
      criterion = criterion, data = x)
}

#' Result of `automarg`
#'
#' @slot best the selected [margfit-class] model.
#' @slot table data frame ranking the families.
#' @slot fits list of the fitted [margfit-class] models of the families that
#'   could be fitted.
#' @slot criterion `"AIC"` or `"BIC"`.
#' @slot data the data used.
#'
#' @exportClass automarg
setClass("automarg", slots = c(best = "margfit", table = "data.frame", fits = "list",
                               criterion = "character", data = "numeric"))

#' @rdname tscop-methods
#' @export
setMethod("show", "automarg", function(object) {
  cat("object class: automarg   (criterion: ", object@criterion, ")\n", sep = "")
  tb <- object@table
  out <- data.frame(family = tb$family, npar = tb$npar, loglik = round(tb$loglik, 2),
                    AIC = round(tb$AIC, 2), BIC = round(tb$BIC, 2), delta = round(tb$delta, 2),
                    stringsAsFactors = FALSE)
  print(out, row.names = FALSE)
  bad <- tb$note != ""
  if (any(bad)) cat("not fitted: ", paste0(tb$family[bad], " (", tb$note[bad], ")", collapse = "; "), "\n", sep = "")
  cat("\nselected model:\n")
  show(object@best)
  invisible(object)
})

#' @rdname tscop-methods
#' @export
setMethod("coef", "automarg", function(object) coef(object@best))
#' @rdname tscop-methods
#' @export
setMethod("logLik", "automarg", function(object) logLik(object@best))

#' @rdname tssim
#' @export
setMethod("tssim", "automarg", function(object, n, ...) {
  if (missing(n)) n <- length(object@data)
  tssim(object@best, n)
})

#' @rdname resid-method
#' @export
setMethod("resid", "automarg", function(object, data, scale = "uniform", nodes = 64, ...) {
  if (missing(data)) resid(object@best, scale = scale) else resid(object@best, data, scale = scale)
})
#' @rdname resid-method
#' @export
setMethod("residuals", "automarg", function(object, data, scale = "uniform", nodes = 64, ...) {
  if (missing(data)) resid(object@best, scale = scale) else resid(object@best, data, scale = scale)
})

#' @rdname marg-functions
#' @export
setMethod("mdens", "automarg", function(object, x, log = FALSE, ...) mdens(object@best, x, log))
#' @rdname marg-functions
#' @export
setMethod("mdf", "automarg", function(object, x, lower.tail = TRUE, ...) mdf(object@best, x, lower.tail))
#' @rdname marg-functions
#' @export
setMethod("mqf", "automarg", function(object, p, lower.tail = TRUE, ...) mqf(object@best, p, lower.tail))
