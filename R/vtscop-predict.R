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

vt_qf_ref <- function(vt, p, Fy, Qy, nodes) {
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

# ---- fast predictive df and qf for ONE conditional law ----------------------
#
# For a nonlinear v-transform, F_U(x) = G(p0(x)) + 1{x > d} p0(x), with
# G(p) = int_p^1 Delta(Q_Y(s)) ds. G depends on the day's conditional law of Y but not on x,
# so it is tabulated once (Delta(Q_Y) at the knots of a Chebyshev-spaced grid, integrated
# through its cubic spline interpolant, then interpolated by cubic Hermite with G' = -Delta(Q_Y)),
# and every later evaluation of the df, and each step of the quantile search, is a lookup plus one
# evaluation of the base df. The quadrature versions (vt_df(), vt_qf_ref()) are the reference.
# Linear v-transforms are closed form, and the quantile function is explicit.

vt_table <- function(vt, Qy, n = 128L) {
  k <- 0:n
  x <- (1 - cos(pi * k / n)) / 2
  y <- pmin(pmax(Qy(pmin(pmax(x, 1e-12), 1 - 1e-12), rep(1L, n + 1L)), 1e-10), 1 - 1e-10)
  h <- vt_delta(vt, y)
  sp <- stats::splinefun(x, h, method = "fmm")
  gl <- gl_nodes(4L)                                  # 4-point Gauss-Legendre on each cell
  a <- x[-(n + 1L)]; w <- diff(x)
  cell <- vapply(seq_along(a), function(j) sum(sp(a[j] + w[j] * gl$s) * gl$w) * w[j], numeric(1))
  G <- c(rev(cumsum(rev(cell))), 0)                   # int_{x_k}^1 of the interpolant
  stats::splinefunH(x, G, -h)
}

vt_df1 <- function(vt, x, Fy, G) {
  d <- vt_fulcrum(vt)
  p0 <- pmin(pmax(Fy(vt_transform_data(vt, x), seq_along(x)), 0), 1)
  I <- if (vt_is_linear(vt)) d * (1 - p0) else G(p0)
  I + (x > d) * p0
}

vt_qf <- function(vt, p, Fy, Qy, nodes) {
  d <- vt_fulcrum(vt)
  if (vt_is_linear(vt)) {                             # explicit inverse
    lower <- p <= d
    y <- Qy(ifelse(lower, 1 - p / d, (p - d) / (1 - d)), rep(1L, length(p)))
    return(ifelse(lower, d * (1 - y), d + (1 - d) * y))
  }
  G <- vt_table(vt, Qy)
  lo <- rep(0, length(p)); hi <- rep(1, length(p))
  for (it in 1:40) {
    mid <- (lo + hi) / 2
    up <- vt_df1(vt, mid, Fy, G) < p
    lo <- ifelse(up, mid, lo); hi <- ifelse(up, hi, mid)
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
  u <- if (is(model@base, "sdvinecop")) check_complete_copula_data(data) else check_copula_data(data, min_n = 1L)
  v <- rep(NA_real_, length(u))
  ok <- !is.na(u)
  v[ok] <- vt_transform_data(model@vt, u[ok])
  vt0 <- model@vt
  pd <- vt_predictive_dist(model@base, v, h)
  switch(type,
    df = vt_df1(vt0, x, pd$Fy, if (vt_is_linear(vt0)) NULL else vt_table(vt0, pd$Qy)),
    qf = vt_qf(vt0, x, pd$Fy, pd$Qy, nodes),
    dens = vt_dens(vt0, x, pd$fy)
  )
}

vt_resid <- function(model, data, scale, nodes) {
  scale <- match.arg(scale, c("uniform", "normal"))
  if (!vt_base_specified(model)) stop("the model has unspecified coefficients; fit it to data first")
  b <- model@base
  vt <- model@vt
  u <- if (is(b, "sdvinecop")) check_complete_copula_data(data) else check_copula_data(data, min_n = 1L)
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

#' @rdname resid-method
#' @export
setMethod("resid", "vtscop", function(object, data, scale = "uniform", nodes = 64, ...) {
  vt_resid(object, data, scale, nodes)
})

#' @rdname resid-method
#' @export
setMethod("residuals", "vtscop", function(object, data, scale = "uniform", nodes = 64, ...) {
  vt_resid(object, data, scale, nodes)
})

vt_fit_data <- function(object) if (is.null(object@template)) object@data else object@template

#' @rdname resid-method
#' @export
setMethod("resid", "vtscopfit", function(object, data, scale = "uniform", nodes = 64, ...) {
  if (missing(data)) data <- vt_fit_data(object)
  vt_resid(object@model, data, scale, nodes)
})

#' @rdname resid-method
#' @export
setMethod("residuals", "vtscopfit", function(object, data, scale = "uniform", nodes = 64, ...) {
  if (missing(data)) data <- vt_fit_data(object)
  vt_resid(object@model, data, scale, nodes)
})

#' @rdname predict-method
#' @export
setMethod("predict", "vtscop", function(object, data, h = 1, type = "df", x = NULL,
                                        nodes = 64, ...) {
  vt_predict(object, data, h, type, x, nodes)
})

#' @rdname predict-method
#' @export
setMethod("predict", "vtscopfit", function(object, data, h = 1, type = "df", x = NULL,
                                           nodes = 64, ...) {
  if (missing(data)) data <- object@data
  vt_predict(object@model, data, h, type, x, nodes)
})
