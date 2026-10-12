#' Empirical (nonparametric) margins
#'
#' A marginal distribution estimated from the data without a parametric family,
#' for the semiparametric approach in which the margin is estimated by the
#' (scaled) empirical distribution function and the copula process is fitted to
#' the resulting pseudo-observations (Bladt, Dias, Han and McNeil, 2026). Two
#' types are available:
#'
#' * `"edf"`: the distribution function is the scaled empirical distribution
#'   function \eqn{\hat F(x) = (1/(n+1)) \sum_i 1\{x_i \le x\}}, which stays
#'   inside (0, 1) (for `proper = FALSE`; below the smallest observation it
#'   equals \eqn{0.5/(n+1)} and above the largest \eqn{(n+0.5)/(n+1)}, an
#'   improper distribution function that is bounded away from 0 and 1;
#'   `proper = TRUE` gives the ordinary empirical distribution function); the
#'   quantile function is the empirical quantile of type `qtype` (see
#'   [stats::quantile()]); the density is a Gaussian kernel density estimate.
#' * `"kernel"`: the fully smoothed margin: the distribution function, density and
#'   quantile function of the Gaussian kernel density estimate.
#'
#' @slot data the finite observations in their original order.
#' @slot sorted the sorted observations.
#' @slot type `"edf"` or `"kernel"`.
#' @slot bw the kernel bandwidth.
#' @slot qtype the type of empirical quantile (1 to 9).
#' @slot proper logical; use the ordinary empirical distribution function
#'   (`type = "edf"` only).
#' @slot template the original input if it was a `ts`, `zoo` or `xts` object,
#'   otherwise `NULL`.
#'
#' @references Bladt, M., Dias, A., Han, J. and McNeil, A. J. (2026).
#'   Semiparametric and parametric distributional forecasting of univariate time
#'   series using non-Gaussian ARMA models based on D-vines. *Canadian Journal
#'   of Statistics*, 54(1), e70033. \doi{10.1002/cjs.70033}
#'
#' @exportClass emarg
setClass("emarg", slots = c(data = "numeric", sorted = "numeric", type = "character", bw = "numeric", qtype = "numeric",
                            proper = "logical", template = "ANY"))

#' Constructor for empirical margins
#'
#' @param x numeric data (non-finite values are dropped).
#' @param type `"edf"` (default) or `"kernel"`; see [emarg-class].
#' @param bw the kernel bandwidth: a number, or the name of a bandwidth selector
#'   (`"SJ"`, the Sheather-Jones method, the default; `"nrd0"`, `"nrd"`, `"ucv"`,
#'   `"bcv"`; see [stats::bw.SJ()]).
#' @param qtype type of the empirical quantile, 1 to 9 (default 7); see
#'   [stats::quantile()].
#' @param proper logical; for `type = "edf"`, use the ordinary empirical
#'   distribution function instead of the scaled one (default `FALSE`).
#'
#' @return An object of class [emarg-class].
#' @export
#' @examples
#' set.seed(1)
#' x <- rt(300, 5)
#' e <- emarg(x)
#' mdf(e, c(-1, 0, 1))
#' mqf(e, c(0.05, 0.5, 0.95))
#' head(resid(e))          # pseudo-observations rank / (n + 1)
#' mdens(emarg(x, "kernel"), 0)
emarg <- function(x, type = c("edf", "kernel"), bw = "SJ", qtype = 7, proper = FALSE) {
  type <- match.arg(type)
  template <- if (inherits(x, c("ts", "zoo", "xts"))) x else NULL
  if (!is.numeric(x)) stop("the data must be numeric", call. = FALSE)
  xs <- as.numeric(x)
  xs <- xs[is.finite(xs)]
  if (length(xs) < 5L) stop("at least 5 finite observations are needed", call. = FALSE)
  if (!(is.numeric(qtype) && length(qtype) == 1L && qtype %in% 1:9)) stop("'qtype' must be an integer from 1 to 9", call. = FALSE)
  if (is.character(bw)) {
    bw <- match.arg(bw, c("SJ", "nrd0", "nrd", "ucv", "bcv"))
    h <- switch(bw, SJ = tryCatch(stats::bw.SJ(xs), error = function(e) stats::bw.nrd0(xs)),
                nrd0 = stats::bw.nrd0(xs), nrd = stats::bw.nrd(xs), ucv = stats::bw.ucv(xs), bcv = stats::bw.bcv(xs))
  } else if (is.numeric(bw) && length(bw) == 1L && bw > 0) {
    h <- bw
  } else {
    stop("'bw' must be a positive number or the name of a bandwidth selector", call. = FALSE)
  }
  new("emarg", data = xs, sorted = sort(xs), type = type, bw = h, qtype = qtype, proper = isTRUE(proper), template = template)
}

# Kernel quantities by chunks of evaluation points (n x m matrices of at most
# about 2e6 entries).
em_kernel <- function(e, x, what) {
  n <- length(e@sorted)
  h <- e@bw
  m <- length(x)
  out <- numeric(m)
  chunk <- max(1L, floor(2e6 / n))
  for (s in seq(1L, m, by = chunk)) {
    i <- s:min(m, s + chunk - 1L)
    u <- outer(x[i], e@sorted, "-") / h
    out[i] <- switch(what,
      cdf = rowMeans(stats::pnorm(u)),
      sf = rowMeans(stats::pnorm(u, lower.tail = FALSE)),
      dens = rowMeans(stats::dnorm(u)) / h
    )
  }
  out
}

em_edf <- function(e, x) {
  n <- length(e@sorted)
  cnt <- findInterval(x, e@sorted) # number of observations <= x
  if (e@proper) return(cnt / n)
  u <- cnt / (n + 1)
  u[x < e@sorted[1L]] <- 0.5 / (n + 1)
  u[x > e@sorted[n]] <- (n + 0.5) / (n + 1)
  u
}

#' @rdname marg-functions
#' @export
setMethod("mdens", "emarg", function(object, x, log = FALSE, ...) {
  d <- em_kernel(object, as.numeric(x), "dens")
  if (log) base::log(d) else d
})

#' @rdname marg-functions
#' @export
setMethod("mdf", "emarg", function(object, x, lower.tail = TRUE, ...) {
  x <- as.numeric(x)
  if (object@type == "kernel") {
    em_kernel(object, x, if (lower.tail) "cdf" else "sf")
  } else {
    u <- em_edf(object, x)
    if (lower.tail) u else 1 - u
  }
})

#' @rdname marg-functions
#' @export
setMethod("mqf", "emarg", function(object, p, lower.tail = TRUE, ...) {
  p <- as.numeric(p)
  if (any(p < 0 | p > 1, na.rm = TRUE)) stop("'p' must be in [0, 1]", call. = FALSE)
  if (!lower.tail) p <- 1 - p
  q0 <- as.numeric(stats::quantile(object@sorted, p, type = object@qtype, na.rm = FALSE))
  if (object@type == "edf") return(q0)
  # kernel: invert the smoothed distribution function by safeguarded Newton steps
  lo <- object@sorted[1L] - 12 * object@bw
  hi <- object@sorted[length(object@sorted)] + 12 * object@bw
  z <- pmin(pmax(q0, lo), hi)
  a <- rep(lo, length(z))
  b <- rep(hi, length(z))
  ok <- !is.na(p) & p > 0 & p < 1
  for (it in 1:60) {
    if (!any(ok)) break
    F <- em_kernel(object, z[ok], "cdf")
    f <- em_kernel(object, z[ok], "dens")
    above <- F > p[ok]
    ai <- a[ok]; bi <- b[ok]
    bi[above] <- z[ok][above]
    ai[!above] <- z[ok][!above]
    a[ok] <- ai; b[ok] <- bi
    zn <- z[ok] - (F - p[ok]) / f
    bad <- !is.finite(zn) | zn <= ai | zn >= bi
    zn[bad] <- 0.5 * (ai[bad] + bi[bad])
    conv <- abs(zn - z[ok]) < 1e-12 * (1 + abs(zn))
    z[ok] <- zn
    ok[which(ok)[conv]] <- FALSE
  }
  z[!is.na(p) & p <= 0] <- -Inf
  z[!is.na(p) & p >= 1] <- Inf
  z
})

#' @describeIn tssim Resampling from an [emarg-class] margin: the bootstrap
#'   (random draws of the observations) for `type = "edf"`, the smoothed
#'   bootstrap (draws plus Gaussian kernel noise) for `type = "kernel"`; `n`
#'   defaults to the number of observations.
#' @export
setMethod("tssim", "emarg", function(object, n, ...) {
  if (missing(n)) n <- length(object@data)
  if (!is.numeric(n) || length(n) != 1L || n < 1 || n != round(n)) stop("'n' must be a positive whole number")
  x <- object@data[sample.int(length(object@data), n, replace = TRUE)]
  if (object@type == "kernel") x <- x + object@bw * stats::rnorm(n)
  x
})

#' @rdname resid-method
#' @export
setMethod("resid", "emarg", function(object, data, scale = "uniform", nodes = 64, ...) {
  scale <- match.arg(scale, c("uniform", "normal"))
  if (missing(data)) data <- if (is.null(object@template)) object@data else object@template
  if (!is.numeric(data)) stop("'data' must be numeric", call. = FALSE)
  u <- mdf(object, as.numeric(data))
  u <- pmin(pmax(u, 1e-12), 1 - 1e-12)
  restore_attrs(if (scale == "uniform") u else stats::qnorm(u), data)
})

#' @rdname resid-method
#' @export
setMethod("residuals", "emarg", function(object, data, scale = "uniform", nodes = 64, ...) {
  if (missing(data)) resid(object, scale = scale) else resid(object, data, scale = scale)
})

#' @rdname tscop-methods
#' @export
setMethod("show", "emarg", function(object) {
  cat("object class: emarg\n")
  cat("type: ", object@type, if (object@type == "edf") paste0(if (object@proper) " (proper)" else " (scaled)"), "\n", sep = "")
  cat("observations: ", length(object@data), "   kernel bandwidth: ", format(object@bw, digits = 4),
      "   quantile type: ", object@qtype, "\n", sep = "")
  invisible(object)
})
