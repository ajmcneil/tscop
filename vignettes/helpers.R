## ---- helpers ----
library(rugarch)

# the fitted d-vine copula of a tscop fit (a vtscop model: its base d-vine)
cop_of <- function(f) {
  m <- if (is(f, "tsmodfit")) f@copfit@model else f@model
  if (is(m, "vtscop")) m@base else m
}
# log-likelihood, number of parameters and AIC of a fitted model, tscop or rugarch
# (rugarch reports the AIC per observation)
ic <- function(f) {
  if (is(f, "uGARCHfit")) {
    c(loglik = likelihood(f), npar = length(coef(f)),
      AIC = infocriteria(f)["Akaike", 1] * length(residuals(f)))
  } else {
    c(loglik = logLik(f)[1], npar = attr(logLik(f), "df"), AIC = AIC(f))
  }
}
# AIC table of a list of fitted copula models, with the effective lag K and the
# number of nonlinear lags L
aic_table <- function(fits) {
  tab <- data.frame(model = names(fits), t(sapply(fits, ic)),
                    K = sapply(fits, function(f) cop_of(f)@K),
                    L = sapply(fits, function(f) cop_of(f)@L), row.names = NULL)
  tab$dAIC <- tab$AIC - min(tab$AIC)
  tab[order(tab$AIC), ]
}

# g describes a GARCH model: the variance model, the innovation distribution and
# the conditional mean model as a string "(p,q)" for the ARMA orders
garch_spec <- function(g) {
  arma <- as.numeric(strsplit(gsub("[()]", "", g$meanmod), ",")[[1]])
  ugarchspec(variance.model = list(model = g$model),
             mean.model = list(armaOrder = arma),
             distribution.model = g$dist)
}
garch_fit <- function(x, g) {
  f <- tryCatch(ugarchfit(garch_spec(g), x, solver = "hybrid"),
                error = function(e) NULL)
  if (is.null(f) || f@fit$convergence != 0) NULL else f
}
# predictive quantiles (one row per observation, one column per level) of the
# observations in `x`, from the coefficients `cf` of a fit; the filter is run over
# the whole series `x`
garch_quantiles <- function(cf, x, g, levels) {
  sp <- garch_spec(g)
  setfixed(sp) <- as.list(cf)
  fl <- ugarchfilter(sp, x)
  s <- as.numeric(sigma(fl))
  mu <- x - as.numeric(residuals(fl))
  skew <- if ("skew" %in% names(cf)) cf[["skew"]] else 1
  shape <- if ("shape" %in% names(cf)) cf[["shape"]] else 5
  mu + outer(s, qdist(g$dist, levels, skew = skew, shape = shape))
}
# rolling re-estimation: refit every H observations and forecast the next block
garch_roll <- function(x, g, start, H = 20, levels = (1:99) / 100, cf0 = NULL) {
  n <- length(x)
  Q <- matrix(NA_real_, n - start + 1, length(levels))
  cf <- cf0
  for (b0 in seq(start, n, by = H)) {
    idx <- b0:min(b0 + H - 1, n)
    if (!is.null(cf0) && b0 > start || is.null(cf)) {
      f <- garch_fit(x[1:(b0 - 1)], g)
      if (!is.null(f)) cf <- coef(f)
    }
    q <- garch_quantiles(cf, x[1:max(idx)], g, levels)
    Q[idx - start + 1, ] <- q[idx, ]
  }
  Q
}
