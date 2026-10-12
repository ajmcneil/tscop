# Prediction and predictive risk measures for tsmod models.
#
# Given the past data x_1, ..., x_t, the predictive distribution of x_{t+h} is
# that of the copula process for u_{t+h} given the past probability integral
# transforms u_s = F(x_s), pushed through the margin:
#   df(x)   = copula-df(F(x)),
#   qf(p)   = F^{-1}(copula-qf(p)),
#   dens(x) = f(x) * copula-dens(F(x)).

tm_predict <- function(model, data, h, type, x, nodes) {
  type <- match.arg(type, c("df", "qf", "dens"))
  if (!is.numeric(data) || anyNA(data) || any(!is.finite(data)) || length(data) < 1L) {
    stop("'data' must be complete finite past observations", call. = FALSE)
  }
  if (is.null(x) || !is.numeric(x) || anyNA(x)) stop("'x' must be numeric values without missing values", call. = FALSE)
  if (type == "qf" && any(x <= 0 | x >= 1)) stop("for type = \"qf\", 'x' must be probabilities strictly inside (0, 1)", call. = FALSE)
  up <- tm_pits(model@margin, data)
  mg <- model@margin
  cop <- model@copula
  switch(type,
    df = predict(cop, up, h = h, type = "df", x = tm_pits(mg, x), nodes = nodes),
    qf = mqf(mg, predict(cop, up, h = h, type = "qf", x = x, nodes = nodes)),
    dens = {
      ux <- tm_pits(mg, x)
      mdens(mg, x) * predict(cop, up, h = h, type = "dens", x = ux, nodes = nodes)
    }
  )
}

#' @rdname predict-method
#' @export
setMethod("predict", "tsmod", function(object, data, h = 1, type = "df", x = NULL, nodes = 64, ...) {
  tm_predict(object, data, h, type, x, nodes)
})

#' @rdname predict-method
#' @export
setMethod("predict", "tsmodfit", function(object, data, h = 1, type = "df", x = NULL, nodes = 64, ...) {
  if (missing(data)) data <- object@data
  tm_predict(object@model, data, h, type, x, nodes)
})

#' Predictive value-at-risk and expected shortfall
#'
#' Risk measures of the predictive distribution of \eqn{x_{t+h}} given the past
#' data of a [tsmod-class] model. For each level `p` and `tail = "lower"` (the
#' left tail, appropriate for returns) the value-at-risk is the predictive
#' quantile \eqn{q_p = F^{-1}_{t+h|t}(p)} and the expected shortfall is
#' \eqn{E[X \mid X \le q_p] = p^{-1}\int_0^p F^{-1}_{t+h|t}(s)\,ds}; for
#' `tail = "upper"` they are \eqn{q_{1-p}} and
#' \eqn{E[X \mid X \ge q_{1-p}] = p^{-1}\int_{1-p}^1 F^{-1}_{t+h|t}(s)\,ds}.
#' To obtain risk measures for losses, apply `tsrisk()` with `tail = "upper"`
#' to a model of the losses (or use the negative returns). The integral is
#' computed by Gauss-Legendre quadrature after the substitution
#' \eqn{s = p e^{-v}}, which handles the singularity of the quantile function in
#' the tail; it requires a finite mean.
#'
#' @param object a [tsmod-class] model with known margin and copula, or a
#'   [tsmodfit-class] object.
#' @param data past observations (complete); for a fitted object the data used in
#'   the fit by default.
#' @param p tail probabilities (default `c(0.01, 0.05)`).
#' @param tail `"lower"` (default) or `"upper"`.
#' @param h horizon (see [predict-method]).
#' @param nodes number of quadrature nodes for the expected shortfall (and for
#'   nonlinear v-transforms in the copula).
#'
#' @return A data frame with columns `p`, `tail`, `h`, `var` (the predictive
#'   quantile) and `es` (the expected shortfall).
#' @export
#' @examples
#' set.seed(1)
#' m <- tsmod(marg("st", c(mu = 0, sigma = 1, df = 5)), arfactcop(lambda = 0.9, w = 0.6))
#' x <- tssim(m, 300)
#' tsrisk(m, x, p = c(0.01, 0.05))
tsrisk <- function(object, data, p = c(0.01, 0.05), tail = c("lower", "upper"), h = 1, nodes = 64) {
  tail <- match.arg(tail)
  if (is(object, "tsmodfit")) {
    if (missing(data)) data <- object@data
    object <- object@model
  }
  if (!is(object, "tsmod")) stop("'object' must be a tsmod or tsmodfit object", call. = FALSE)
  if (!is.numeric(p) || anyNA(p) || any(p <= 0 | p >= 1)) stop("'p' must be probabilities strictly inside (0, 1)", call. = FALSE)
  g <- gl_nodes(nodes)
  out <- data.frame(p = p, tail = tail, h = h, var = NA_real_, es = NA_real_)
  for (i in seq_along(p)) {
    pi_ <- p[i]
    level <- if (tail == "lower") pi_ else 1 - pi_
    out$var[i] <- tm_predict(object, data, h, "qf", level, nodes)
    # integral of the predictive quantile function over the tail: s = p exp(-v), with
    # the range of v limited so that probabilities stay distinguishable from 0 and 1
    vmax <- min(40, log(pi_ / 1e-14))
    v <- g$s * vmax
    w <- g$w * vmax
    s <- pi_ * exp(-v)
    probs <- if (tail == "lower") s else 1 - s
    q <- tm_predict(object, data, h, "qf", probs, nodes)
    out$es[i] <- sum(w * q * exp(-v))
  }
  out
}
