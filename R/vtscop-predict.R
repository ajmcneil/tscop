# Residuals and prediction for vtscop models.
#
# Let Y = V(U) be the base process. Given the past, the conditional law of
# U_{n+h} is obtained from that of Y_{n+h} by stochastic inversion of the
# v-transform. With Delta(y) the probability of the lower branch of
# V^{-1}(y), and p0 = F_Y(V(x)),
#   F_U(x) = int_{p0}^1 Delta(Q_Y(p)) dp + 1{x > fulcrum} p0,
# which is closed form (Delta constant) for linear V and a one-dimensional
# Gauss-Legendre integral on the probability scale otherwise.

.gl_cache <- new.env(parent = emptyenv())

# Gauss-Legendre nodes and weights on (0, 1) (Golub-Welsch).
gl_nodes <- function(n) {
  key <- as.character(n)
  if (!is.null(.gl_cache[[key]])) return(.gl_cache[[key]])
  i <- seq_len(n - 1L)
  b <- i / sqrt(4 * i^2 - 1)
  J <- matrix(0, n, n)
  J[cbind(i, i + 1L)] <- b
  J[cbind(i + 1L, i)] <- b
  e <- eigen(J, symmetric = TRUE)
  o <- order(e$values)
  out <- list(s = (e$values[o] + 1) / 2, w = e$vectors[1L, o]^2)
  assign(key, out, envir = .gl_cache)
  out
}

vt_fulcrum <- function(vt) if ("delta" %in% names(vt@pars)) vt@pars[["delta"]] else 0.5
vt_is_linear <- function(vt) vt@name %in% c("vlinear", "vsymmetric")

# Lower-branch probabilities Delta(y).
vt_delta <- function(vt, y) attr(udp::udpinverse(vt, y, prob = TRUE), "prob")[, 1L]

# Distribution function of U at x. Fy(y, i) and Qy(p, i) are the base's
# conditional df and quantile function; i indexes the element of x (and so
# the conditional law) that each value belongs to.
vt_df <- function(vt, x, Fy, Qy, nodes) {
  nx <- length(x)
  d <- vt_fulcrum(vt)
  y0 <- vt_transform_data(vt, x)
  p0 <- pmin(pmax(Fy(y0, seq_len(nx)), 0), 1)
  if (vt_is_linear(vt)) {
    I <- d * (1 - p0)
  } else {
    g <- gl_nodes(nodes)
    P <- p0 + outer(1 - p0, g$s)
    Y <- pmin(pmax(Qy(as.vector(P), rep(seq_len(nx), times = nodes)), 1e-10), 1 - 1e-10)
    D <- matrix(vt_delta(vt, Y), nx)
    I <- (1 - p0) * as.vector(D %*% g$w)
  }
  I + (x > d) * p0
}

vt_qf <- function(vt, p, Fy, Qy, nodes) {
  lo <- rep(0, length(p))
  hi <- rep(1, length(p))
  for (it in 1:40) {
    mid <- (lo + hi) / 2
    f <- vt_df(vt, mid, Fy, Qy, nodes)
    up <- f < p
    lo <- ifelse(up, mid, lo)
    hi <- ifelse(up, hi, mid)
  }
  (lo + hi) / 2
}

# Density of U at x; fy(y, i) is the base's conditional density.
vt_dens <- function(vt, x, fy) {
  d <- vt_fulcrum(vt)
  y <- vt_transform_data(vt, x)
  pr <- attr(udp::udpinverse(vt, y, prob = TRUE), "prob")
  bp <- ifelse(x <= d, pr[, 1L], pr[, 2L])
  fy(y, seq_along(x)) * abs(udp::udpderiv(vt, x)) * bp
}

# Conditional law of Y_{n+h} given the (clamped) transformed past v.
vt_predictive_dist <- function(base, v, h) {
  if (is(base, "arfactcop")) {
    m <- afc_predict(base, v, h, "moments", NULL)
    mu <- m$mean
    s <- m$sd
    list(
      Fy = function(y, i) stats::pnorm((stats::qnorm(y) - mu[1L]) / s[1L]),
      Qy = function(p, i) stats::pnorm(mu[1L] + s[1L] * stats::qnorm(p)),
      fy = function(y, i) stats::dnorm(stats::qnorm(y), mu[1L], s[1L]) / stats::dnorm(stats::qnorm(y))
    )
  } else {
    if (h != 1) stop("only h = 1 is implemented for a d-vine base")
    pairs <- dv_pairs_or_stop(base)
    list(
      Fy = function(y, i) dv_predictive(pairs, v, y, "df"),
      Qy = function(p, i) dv_predictive(pairs, v, p, "qf"),
      fy = function(y, i) dv_predictive(pairs, v, y, "dens")
    )
  }
}

vt_predict <- function(model, data, h, type, x, nodes) {
  type <- match.arg(type, c("df", "qf", "dens"))
  if (!is.numeric(h) || length(h) != 1L || h < 1 || h != round(h)) stop("'h' must be a positive whole number")
  if (is.null(x) || !is.numeric(x) || anyNA(x) || any(x <= 0 | x >= 1)) {
    stop("'x' must be numeric values strictly inside (0, 1)")
  }
  if (!vt_base_specified(model)) stop("the model has unspecified coefficients; fit it to data first")
  u <- if (is(model@base, "dvinecop")) check_complete_copula_data(data) else check_copula_data(data, min_n = 1L)
  v <- rep(NA_real_, length(u))
  ok <- !is.na(u)
  v[ok] <- vt_transform_data(model@vt, u[ok])
  vt0 <- model@vt
  pd <- vt_predictive_dist(model@base, v, h)
  switch(type,
    df = vt_df(vt0, x, pd$Fy, pd$Qy, nodes),
    qf = vt_qf(vt0, x, pd$Fy, pd$Qy, nodes),
    dens = vt_dens(vt0, x, pd$fy)
  )
}

vt_resid <- function(model, data, scale, nodes) {
  scale <- match.arg(scale, c("uniform", "normal"))
  if (!vt_base_specified(model)) stop("the model has unspecified coefficients; fit it to data first")
  b <- model@base
  vt <- model@vt
  u <- if (is(b, "dvinecop")) check_complete_copula_data(data) else check_copula_data(data, min_n = 1L)
  n <- length(u)
  ok <- !is.na(u)
  v <- rep(NA_real_, n)
  v[ok] <- vt_transform_data(vt, u[ok])
  r <- rep(NA_real_, n)
  linear <- vt_is_linear(vt)
  if (is(b, "arfactcop")) {
    p <- afc_check_known(b)
    f <- afc_filter(stats::qnorm(v), p$lambda, p$w, p$w0)
    s <- sqrt(f$F)
    mu <- stats::qnorm(v) - f$v
    Fy <- function(y, i) stats::pnorm((stats::qnorm(y) - mu[i]) / s[i])
    Qy <- function(q, i) stats::pnorm(mu[i] + s[i] * stats::qnorm(q))
    r[ok] <- vt_df(vt, u[ok], function(y, i) Fy(y, which(ok)[i]),
                   function(q, i) Qy(q, which(ok)[i]), nodes)
  } else {
    pairs <- dv_pairs_or_stop(b)
    if (linear) {
      p0 <- dv_rosenblatt(pairs, v)
      d <- vt_fulcrum(vt)
      r <- d * (1 - p0) + (u > d) * p0
    } else {
      K <- length(pairs)
      for (t in seq_len(n)) {
        past <- if (t > 1L) v[max(1L, t - K):(t - 1L)] else numeric(0)
        r[t] <- vt_df(vt, u[t],
                      function(y, i) dv_predictive(pairs, past, y, "df"),
                      function(q, i) dv_predictive(pairs, past, q, "qf"), nodes)
      }
    }
  }
  r <- pmin(pmax(r, 0), 1)
  restore_attrs(if (scale == "uniform") r else stats::qnorm(r), data)
}

#' Residuals of process-level v-transform copula models
#'
#' The residuals are the conditional probability integral transforms
#' \eqn{F(u_t \mid u_{t-1}, \ldots)} of the series, iid uniform under the
#' model. The conditional distribution of \eqn{U_t} is that of the base
#' process at \eqn{V(U_t)} passed through the stochastic inverse of the
#' v-transform: closed form for linear v-transforms, otherwise a
#' one-dimensional Gauss-Legendre integral. For a d-vine base with a
#' nonlinear v-transform the residuals are computed one time point at a time
#' and can be slow.
#'
#' @param object a [vtscop-class] model with known coefficients or a
#'   [vtscopfit-class] object.
#' @param data copula data in (0, 1); for a fitted object the data used in the
#'   fit by default.
#' @param scale `"uniform"` (default) or `"normal"`.
#' @param nodes number of quadrature nodes for nonlinear v-transforms.
#' @param ... not used.
#'
#' @return A numeric vector the length of the data.
#' @name resid-vtscop
#' @aliases resid,vtscop-method residuals,vtscop-method
NULL

#' @rdname resid-vtscop
#' @export
setMethod("resid", "vtscop", function(object, data, scale = "uniform", nodes = 64, ...) {
  vt_resid(object, data, scale, nodes)
})

#' @rdname resid-vtscop
#' @export
setMethod("residuals", "vtscop", function(object, data, scale = "uniform", nodes = 64, ...) {
  vt_resid(object, data, scale, nodes)
})

vt_fit_data <- function(object) if (is.null(object@template)) object@data else object@template

#' @rdname resid-vtscop
#' @export
setMethod("resid", "vtscopfit", function(object, data, scale = "uniform", nodes = 64, ...) {
  if (missing(data)) data <- vt_fit_data(object)
  vt_resid(object@model, data, scale, nodes)
})

#' @rdname resid-vtscop
#' @export
setMethod("residuals", "vtscopfit", function(object, data, scale = "uniform", nodes = 64, ...) {
  if (missing(data)) data <- vt_fit_data(object)
  vt_resid(object@model, data, scale, nodes)
})

#' Prediction for process-level v-transform copula models
#'
#' Predictive distribution of \eqn{u_{n+h}} given the data: the distribution
#' function (`type = "df"`), quantile function (`type = "qf"`) or density
#' (`type = "dens"`) evaluated at `x`. The base process's conditional law of
#' \eqn{V_{n+h}} given \eqn{V(u_1), \ldots, V(u_n)} is passed through the
#' stochastic inverse of the v-transform. Horizons `h > 1` are available for an
#' `arfactcop` base only; a d-vine base has `h = 1`.
#'
#' @param object a [vtscop-class] model with known coefficients or a
#'   [vtscopfit-class] object.
#' @param data past copula data in (0, 1); for a fitted object the data used in
#'   the fit by default.
#' @param h horizon.
#' @param type `"df"`, `"qf"` or `"dens"`.
#' @param x points in (0, 1) at which to evaluate the function (for `"qf"`,
#'   probabilities).
#' @param nodes number of quadrature nodes for nonlinear v-transforms.
#' @param ... not used.
#'
#' @return A numeric vector the length of `x`.
#' @name predict-vtscop
#' @aliases predict,vtscop-method
NULL

#' @rdname predict-vtscop
#' @export
setMethod("predict", "vtscop", function(object, data, h = 1, type = "df", x = NULL,
                                        nodes = 64, ...) {
  vt_predict(object, data, h, type, x, nodes)
})

#' @rdname predict-vtscop
#' @export
setMethod("predict", "vtscopfit", function(object, data, h = 1, type = "df", x = NULL,
                                           nodes = 64, ...) {
  if (missing(data)) data <- object@data
  vt_predict(object@model, data, h, type, x, nodes)
})
