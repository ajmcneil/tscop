#' AR-factor copula processes
#'
#' The Gaussian copula of a sum of independent AR(1) "factors" plus white noise.
#' The normal scores \eqn{x_t = \Phi^{-1}(u_t)} are
#' \deqn{x_t = \sum_{i=1}^m \sqrt{w_i}\, Y_{i,t} + \sqrt{w_0}\,\varepsilon_t,}
#' where the \eqn{Y_i} are independent stationary unit-variance Gaussian AR(1)
#' processes with coefficients \eqn{\lambda_1 > \dots > \lambda_m}, in (0, 1),
#' \eqn{\varepsilon_t} is white noise, the weights \eqn{w_i \ge 0} and
#' \eqn{w_0 = 1 - \sum w_i \ge 0}, so that \eqn{\mathrm{Var}(x_t) = 1}. Its
#' autocorrelations are \eqn{\rho_k = \sum_i w_i \lambda_i^k} for
#' \eqn{k \ge 1}; the process is a Gaussian ARMA(m, m) process (ARMA(m, m - 1)
#' without noise) with real positive AR eigenvalues. Its partial
#' autocorrelations are non-negative at all lags in every numerical check
#' (random draws with up to four factors; this is a conjecture, not a
#' theorem), which makes it a suitable kpacf for the inverse-v-transformed
#' base copulas of volatility models.
#'
#' The free parameters are the `m` rates and the `m` weights (the noise share
#' is the remainder); with `noise = FALSE` the weights sum to one and there
#' are `2m - 1` free parameters (the coefficient vector then contains the
#' first `m - 1` weights).
#'
#' @slot name short description.
#' @slot m number of factors.
#' @slot noise logical: is there a white-noise share \eqn{w_0}?
#' @slot lambda the AR(1) coefficients, decreasing (`NA` in a structure).
#' @slot w the factor weights (`NA` in a structure).
#'
#' @references Bladt, M. and McNeil, A. J. (2022). Time series with infinite-order partial
#'   copula dependence. *Dependence Modeling*, 10(1), 87-107.
#'   \doi{10.1515/demo-2022-0105}
#'
#'   Bladt, M., Dias, A., Han, J. and McNeil, A. J. (2026). Semiparametric and
#'   parametric distributional forecasting of univariate time series using
#'   non-Gaussian ARMA models based on D-vines. *Canadian Journal of Statistics*,
#'   54(1), e70033. \doi{10.1002/cjs.70033}
#'
#'   McNeil, A. J. (2021). Modelling volatile time series with v-transforms and
#'   copulas. *Risks*, 9(1), 14. \doi{10.3390/risks9010014}
#'
#' @exportClass arfactcop
setClass("arfactcop", contains = "gausscop", slots = c(
  name = "character", m = "integer", noise = "logical", lambda = "numeric", w = "numeric"
))

setValidity("arfactcop", function(object) {
  m <- object@m
  if (length(m) != 1L || is.na(m) || m < 1L) return("'m' must be a positive integer")
  if (length(object@lambda) != m || length(object@w) != m) {
    return("'lambda' and 'w' must have length m")
  }
  if (length(object@noise) != 1L || is.na(object@noise)) return("'noise' must be TRUE or FALSE")
  if (!anyNA(c(object@lambda, object@w))) {
    if (any(object@lambda <= 0 | object@lambda >= 1)) return("the AR(1) coefficients must lie in (0, 1)")
    if (m > 1L && any(diff(object@lambda) >= 0)) {
      return("the AR(1) coefficients must be strictly decreasing (distinct)")
    }
    if (any(object@w < 0)) return("the weights must be non-negative")
    if (sum(object@w) > 1 + 1e-10) return("the weights must sum to at most one")
    if (!object@noise && abs(sum(object@w) - 1) > 1e-8) {
      return("without noise the weights must sum to one")
    }
  }
  TRUE
})

#' Constructor for AR-factor copula processes
#'
#' Creates an [arfactcop-class] object, either as a known model by giving the
#' coefficients, or as a structure to be fitted by giving the number of
#' factors.
#'
#' @param lambda the AR(1) coefficients of the factors, in (0, 1) (any order;
#'   they are sorted in decreasing order together with their weights).
#' @param w the weights of the factors, non-negative, with sum at most one
#'   (equal to one if `noise = FALSE`).
#' @param order the number of factors `m` of a structure to be fitted.
#' @param noise logical: include a white-noise component with weight
#'   `1 - sum(w)`? Default `TRUE`.
#'
#' @return An object of class [arfactcop-class].
#' @references Bladt, M. and McNeil, A. J. (2022). Time series with infinite-order partial
#'   copula dependence. *Dependence Modeling*, 10(1), 87-107.
#'   \doi{10.1515/demo-2022-0105}
#'
#'   Bladt, M., Dias, A., Han, J. and McNeil, A. J. (2026). Semiparametric and
#'   parametric distributional forecasting of univariate time series using
#'   non-Gaussian ARMA models based on D-vines. *Canadian Journal of Statistics*,
#'   54(1), e70033. \doi{10.1002/cjs.70033}
#'
#'   McNeil, A. J. (2021). Modelling volatile time series with v-transforms and
#'   copulas. *Risks*, 9(1), 14. \doi{10.3390/risks9010014}
#'
#' @export
#'
#' @examples
#' arfactcop(lambda = c(0.98, 0.7), w = c(0.3, 0.4))   # noise share 0.3
#' arfactcop(lambda = 0.9, w = 1, noise = FALSE)       # an AR(1)
#' arfactcop(order = 2)                                # to be fitted
arfactcop <- function(lambda = numeric(0), w = numeric(0), order = NULL, noise = TRUE) {
  if (!is.logical(noise) || length(noise) != 1L || is.na(noise)) stop("'noise' must be TRUE or FALSE")
  if (!is.null(order)) {
    if (length(lambda) > 0L || length(w) > 0L) stop("give either coefficients or 'order', not both")
    if (!is.numeric(order) || length(order) != 1L || is.na(order) || order < 1 || order != round(order)) {
      stop("'order' must be a positive whole number")
    }
    m <- as.integer(order)
    lambda <- rep(NA_real_, m)
    w <- rep(NA_real_, m)
  } else {
    if (length(lambda) < 1L) stop("specify 'lambda' and 'w', or 'order'")
    if (!is.numeric(lambda) || !is.numeric(w) || anyNA(lambda) || anyNA(w)) {
      stop("'lambda' and 'w' must be numeric without missing values")
    }
    if (length(lambda) != length(w)) stop("'lambda' and 'w' must have the same length")
    m <- length(lambda)
    o <- order(lambda, decreasing = TRUE)
    lambda <- as.numeric(lambda[o])
    w <- as.numeric(w[o])
    if (!noise && abs(sum(w) - 1) <= 1e-8) w <- w / sum(w)
  }
  new("arfactcop",
    name = sprintf("ARFACT(%d)%s", m, if (noise) "" else " without noise"),
    m = m, noise = noise, lambda = lambda, w = w
  )
}

# Natural (free) coefficients: rates, then the weights (the first m - 1 without
# noise).
afc_coef <- function(object) {
  m <- object@m
  w <- if (object@noise) object@w else object@w[seq_len(m - 1L)]
  out <- c(object@lambda, w)
  names(out) <- c(paste0("lambda", seq_len(m)), if (length(w) > 0L) paste0("w", seq_along(w)))
  out
}

#' @rdname tscop-methods
#' @export
setMethod("coef", "arfactcop", function(object) afc_coef(object))

#' @rdname tscop-methods
#' @export
setMethod("show", "arfactcop", function(object) {
  cat("object class: arfactcop\n")
  cat("name: ", object@name, "\n", sep = "")
  if (anyNA(c(object@lambda, object@w))) {
    cat("parameters (not yet estimated):\n")
  } else {
    cat("parameters:\n")
  }
  print(afc_coef(object))
  if (!anyNA(c(object@lambda, object@w)) && object@noise) {
    cat("noise share w0: ", format(1 - sum(object@w), digits = 4), "\n", sep = "")
  }
  invisible(object)
})

# ---- autocorrelations and partial autocorrelations ------------------------

# Autocorrelations at lags 1..K of the factor model.
afc_acf <- function(lambda, w, K) {
  if (K < 1L) return(numeric(0))
  colSums(w * outer(lambda, seq_len(K), `^`))
}

# Partial autocorrelations from autocorrelations rho_1, ..., rho_L by the
# Durbin-Levinson recursion.
acf_to_pacf <- function(rho) {
  L <- length(rho)
  pacf <- numeric(L)
  phi <- numeric(0)
  for (k in seq_len(L)) {
    pk <- if (k == 1L) rho[1L] else
      (rho[k] - sum(phi * rho[(k - 1L):1L])) / (1 - sum(phi * rho[1:(k - 1L)]))
    pacf[k] <- pk
    phi <- c(phi - pk * rev(phi), pk)
  }
  pacf
}

# PACF of the factor model; rounding-level negative values (the PACF is
# non-negative) are set to zero.
afc_pacf <- function(lambda, w, K) {
  if (K < 1L) return(numeric(0))
  a <- acf_to_pacf(afc_acf(lambda, w, K))
  a[a < 0 & a > -1e-10] <- 0
  a
}

# ---- free-scale parametrisation ------------------------------------------

afc_nfree <- function(m, noise) if (noise) 2L * m else 2L * m - 1L

# Unconstrained vector -> list(lambda, w, w0). Rates: lambda_1 =
# logistic(z_1), lambda_i = lambda_{i-1} logistic(z_i), so 1 > lambda_1 > ...
# > lambda_m > 0. Weights: softmax over (u_1, ..., u_m, 0) with the noise as
# the reference category, or over (u_1, ..., u_{m-1}, 0) without noise.
afc_free_to_par <- function(z, m, noise) {
  lam <- numeric(m)
  lam[1L] <- stats::plogis(z[1L])
  if (m > 1L) for (i in 2:m) lam[i] <- lam[i - 1L] * stats::plogis(z[i])
  u <- if (length(z) > m) z[(m + 1L):length(z)] else numeric(0)
  e <- c(u, 0)
  e <- exp(e - max(e))
  p <- e / sum(e)
  if (noise) {
    list(lambda = lam, w = p[seq_len(m)], w0 = p[m + 1L])
  } else {
    list(lambda = lam, w = p, w0 = 0)
  }
}

afc_par_to_free <- function(lambda, w, noise) {
  m <- length(lambda)
  z <- numeric(m)
  z[1L] <- stats::qlogis(lambda[1L])
  if (m > 1L) for (i in 2:m) z[i] <- stats::qlogis(lambda[i] / lambda[i - 1L])
  tiny <- 1e-8
  if (noise) {
    w0 <- max(1 - sum(w), tiny)
    u <- log(pmax(w, tiny) / w0)
  } else {
    u <- if (m > 1L) log(pmax(w[-m], tiny) / max(w[m], tiny)) else numeric(0)
  }
  c(z, u)
}

# ---- the gausscop interface -------------------------------------------------

setMethod("gc_is_specified", "arfactcop", function(G) !anyNA(c(G@lambda, G@w)))

setMethod("gc_describe", "arfactcop", function(G) {
  paste0(G@m, "-factor AR(1) model", if (!G@noise) " without noise" else "")
})

setMethod("gc_pacf", "arfactcop", function(G, K) {
  if (!gc_is_specified(G)) stop("the model has unspecified coefficients; fit it to data first")
  afc_pacf(G@lambda, G@w, K)
})

setMethod("gc_effective_lag", "arfactcop", function(G, K0, tautol, maxlag) {
  if (G@m == 1L && !G@noise) return(dv_finish_lag(1L, "ar", K0, maxlag))
  dv_effective_lag_blocks(function(L) afc_pacf(G@lambda, G@w, L), K0, tautol, maxlag)
})

setMethod("gc_structure", "arfactcop", function(G) list(class = "arfactcop", m = G@m, noise = G@noise))

setMethod("gc_from_coef", "arfactcop", function(G, v) {
  m <- G@m
  v <- as.numeric(v)
  lam <- v[seq_len(m)]
  w <- if (G@noise) v[m + seq_len(m)] else {
    ww <- if (m > 1L) v[m + seq_len(m - 1L)] else numeric(0)
    c(ww, 1 - sum(ww))
  }
  arfactcop(lambda = lam, w = w, noise = G@noise)
})

setMethod("gc_free_ops", "arfactcop", function(G) {
  m <- G@m
  noise <- G@noise
  list(
    nfree = afc_nfree(m, noise),
    zbound = 8,
    to_free = function(G2) afc_par_to_free(G2@lambda, G2@w, noise),
    pacf = function(z, K) {
      p <- afc_free_to_par(z, m, noise)
      afc_pacf(p$lambda, p$w, K)
    },
    effective_lag = function(z, K0, tautol, maxlag) {
      p <- afc_free_to_par(z, m, noise)
      if (m == 1L && !noise) return(dv_finish_lag(1L, "ar", K0, maxlag))
      dv_effective_lag_blocks(function(L) afc_pacf(p$lambda, p$w, L), K0, tautol, maxlag)
    },
    coef = function(z) {
      p <- afc_free_to_par(z, m, noise)
      c(p$lambda, if (noise) p$w else p$w[seq_len(m - 1L)])
    },
    build = function(z) {
      p <- afc_free_to_par(z, m, noise)
      arfactcop(lambda = p$lambda, w = p$w, noise = noise)
    }
  )
})

# ---- exact Gaussian likelihood by the Kalman filter ------------------------

# Kalman filter for the normal scores x (NA allowed) of the factor model: the
# state is the vector of the m independent unit-variance AR(1) components,
# x_t = sum_i sqrt(w_i) a_i,t + sqrt(w0) e_t. Returns the predictive means
# 'mu' and variances 'F', the innovations 'v', the Gaussian log-likelihood and
# the predictive state 'a' and covariance 'P' for time n + 1. The variance of
# x_t is exactly one, so no scale parameter is needed.
afc_filter <- function(x, lambda, w, w0, tol = 1e-12) {
  n <- length(x)
  m <- length(lambda)
  zv <- sqrt(w)
  ll2 <- outer(lambda, lambda)
  Q <- diag(1 - lambda^2, m)
  a <- numeric(m)
  P <- diag(m)
  mu <- numeric(n)
  Fv <- numeric(n)
  v <- rep(NA_real_, n)
  ll <- 0
  steady <- FALSE
  Kg <- NULL
  for (t in seq_len(n)) {
    Pz <- as.numeric(P %*% zv)
    Ft <- sum(zv * Pz) + w0
    mu[t] <- sum(zv * a)
    Fv[t] <- Ft
    if (is.na(x[t])) {
      a <- lambda * a
      P <- ll2 * P + Q
      steady <- FALSE
    } else {
      vt <- x[t] - mu[t]
      v[t] <- vt
      ll <- ll - 0.5 * (log(2 * pi) + log(Ft) + vt^2 / Ft)
      if (!steady) Kg <- Pz / Ft
      a <- lambda * (a + Kg * vt)
      if (!steady) {
        Pn <- ll2 * (P - tcrossprod(Pz) / Ft) + Q
        steady <- max(abs(Pn - P)) < tol
        P <- Pn
      }
    }
  }
  list(mu = mu, F = Fv, v = v, loglik = ll, a = a, P = P)
}

# Gaussian log-likelihood of the normal scores: the compiled filter
# (src/afc_kalman.cpp) unless options(tscop.cpp = FALSE), else the R filter.
afc_loglik <- function(x, lambda, w, w0, tol = 1e-12) {
  if (isTRUE(getOption("tscop.cpp", TRUE))) afc_loglik_cpp(as.numeric(x), lambda, w, w0, tol)
  else afc_filter(x, lambda, w, w0, tol)$loglik
}

# Jacobian of the natural parameters (lambda, w, w0) with respect to the free
# vector z, by central differences (the map is smooth and cheap).
afc_free_jacobian <- function(z, m, noise) {
  nz <- length(z)
  Jl <- Jw <- matrix(0, m, nz)
  Jw0 <- numeric(nz)
  for (k in seq_len(nz)) {
    h <- 1e-6 * max(1, abs(z[k]))
    zp <- z; zp[k] <- zp[k] + h
    zm <- z; zm[k] <- zm[k] - h
    pp <- afc_free_to_par(zp, m, noise)
    pm <- afc_free_to_par(zm, m, noise)
    Jl[, k] <- (pp$lambda - pm$lambda) / (2 * h)
    Jw[, k] <- (pp$w - pm$w) / (2 * h)
    Jw0[k] <- (pp$w0 - pm$w0) / (2 * h)
  }
  list(Jl = Jl, Jw = Jw, Jw0 = Jw0)
}

# Gradient of afc_negll (the negative Gaussian log-likelihood of the normal scores
# x) from the compiled filter. The directions are `before` extra data directions,
# then the free parameters z, then `after` extra data directions; dx is the n x
# (before + length(z) + after) matrix of the derivatives of x along them (NULL for
# constant data). NULL if the gradient is unavailable (options, bounds, more than
# 12 directions).
afc_negll_grad <- function(z, m, noise, x, dx = NULL, before = 0L, after = 0L) {
  if (!isTRUE(getOption("tscop.cpp", TRUE)) || !isTRUE(getOption("tscop.cppgrad", TRUE))) return(NULL)
  nz <- length(z)
  P <- before + nz + after
  if (nz < 1L || P > 12L || any(!is.finite(z)) || any(abs(z) > 8 + 1e-8)) return(NULL)
  p <- afc_free_to_par(z, m, noise)
  J <- afc_free_jacobian(z, m, noise)
  pad <- function(M) cbind(matrix(0, nrow(M), before), M, matrix(0, nrow(M), after))
  r <- afc_loglik_grad_cpp(as.numeric(x), if (is.null(dx)) matrix(0, 0, P) else dx, p$lambda, pad(J$Jl), p$w, pad(J$Jw), p$w0,
                           c(numeric(before), J$Jw0, numeric(after)), 1e-12)
  if (!is.finite(r$value) || any(!is.finite(r$gradient))) return(NULL)
  -r$gradient
}

# Gradient for optim(): exact where available, otherwise central differences of afc_negll.
afc_negll_gr <- function(z, m, noise, x) {
  g <- afc_negll_grad(z, m, noise, x)
  if (!is.null(g)) return(g)
  vapply(seq_along(z), function(k) {
    h <- 1e-6
    zp <- z; zp[k] <- zp[k] + h
    zm <- z; zm[k] <- zm[k] - h
    (afc_negll(zp, m, noise, x) - afc_negll(zm, m, noise, x)) / (2 * h)
  }, numeric(1))
}

# Copula log-likelihood of u under a factor model.
afc_copula_loglik <- function(u, lambda, w, w0) {
  x <- stats::qnorm(u)
  afc_filter(x, lambda, w, w0)$loglik - sum(stats::dnorm(x, log = TRUE), na.rm = TRUE)
}

# Negative Gaussian log-likelihood on the free scale; a large value outside the
# bounds or on numerical failure.
afc_negll <- function(z, m, noise, x) {
  if (any(!is.finite(z)) || any(abs(z) > 8 + 1e-8)) return(1e10)   # (tolerance: a fit on the bound rebuilds to 8 + 1e-13)
  p <- afc_free_to_par(z, m, noise)
  val <- tryCatch(-afc_loglik(x, p$lambda, p$w, p$w0), error = function(e) NA_real_)
  if (is.finite(val)) val else 1e10
}

# Default starting values (natural parameters) for m factors.
afc_default_start <- function(m, noise) {
  lam <- if (m == 1L) 0.9 else 0.98 - 0.5 * (seq_len(m) - 1L) / m
  w <- rep(if (noise) 0.6 / m else 1 / m, m)
  afc_par_to_free(lam, w, noise)
}

# Least-squares fit of the factor autocorrelations to the sample
# autocorrelations of the normal scores; a starting vector on the free scale.
afc_ls_start <- function(x, m, noise) {
  L <- min(25L, max(5L, floor(sum(!is.na(x)) / 5)))
  rho_hat <- as.numeric(stats::acf(x, lag.max = L, plot = FALSE, na.action = stats::na.pass)$acf)[-1L]
  obj <- function(z) {
    if (any(abs(z) > 8)) return(1e6)
    p <- afc_free_to_par(z, m, noise)
    sum((rho_hat - afc_acf(p$lambda, p$w, L))^2)
  }
  z0 <- afc_default_start(m, noise)
  o <- stats::optim(z0, obj, method = "Nelder-Mead",
                    control = list(maxit = 2000, reltol = 1e-10, warn.1d.NelderMead = FALSE))
  o$par
}

# Maximum likelihood of the Gaussian factor copula for copula data u.
arfactcop_mle <- function(u, object, init, control) {
  x <- stats::qnorm(u)
  m <- object@m
  noise <- object@noise
  starts <- list()
  if (!is.null(init)) {
    G <- if (is(init, "arfactcop")) init else gc_from_coef(object, init)
    starts[[1L]] <- afc_par_to_free(G@lambda, G@w, noise)
  } else if (gc_is_specified(object)) {
    starts[[1L]] <- afc_par_to_free(object@lambda, object@w, noise)
  } else {
    starts[[1L]] <- afc_ls_start(x, m, noise)
    starts[[2L]] <- afc_default_start(m, noise)
  }
  vals <- vapply(starts, afc_negll, numeric(1), m = m, noise = noise, x = x)
  if (all(vals >= 1e10)) stop("the starting values are infeasible for these data")
  z0 <- starts[[which.min(vals)]]
  ctrl <- modifyList(list(maxit = 500, reltol = 1e-10, ndeps = rep(1e-6, length(z0))), control)
  opt <- stats::optim(z0, afc_negll, gr = afc_negll_gr, m = m, noise = noise, x = x, method = "BFGS", control = ctrl)
  if (opt$convergence != 0L || opt$value >= 1e10) {
    alt <- stats::optim(opt$par, afc_negll, m = m, noise = noise, x = x, method = "Nelder-Mead",
                        control = list(maxit = 5000, reltol = 1e-10, warn.1d.NelderMead = FALSE))
    if (alt$value < opt$value) opt <- alt
  }
  p <- afc_free_to_par(opt$par, m, noise)
  model <- arfactcop(lambda = p$lambda, w = p$w, noise = noise)
  start <- tryCatch({
    ps <- afc_free_to_par(z0, m, noise)
    arfactcop(lambda = ps$lambda, w = ps$w, noise = noise)
  }, error = function(e) model)
  list(model = model, start = start, opt = opt, x = x, z = opt$par)
}

setMethod("gc_prefit", "arfactcop", function(G, U, control) arfactcop_mle(U, G, NULL, control)$model)

# ---- conversion to an equivalent ARMA copula ---------------------------------

#' Convert an AR-factor copula to the equivalent ARMA copula
#'
#' The Gaussian copula of independent AR(1) factors plus noise has the
#' autocorrelations of a Gaussian ARMA(m, m) process (ARMA(m, m - 1) without
#' noise): the AR polynomial is \eqn{\prod_i (1 - \lambda_i B)} and the MA
#' polynomial follows from the spectral factorisation (Fejér–Riesz) of the
#' sum of the factors' spectral densities, choosing the invertible root. This
#' function returns that [armacop-class] model, so the implied ARMA
#' coefficients can be inspected, or the model passed to ARMA software. The two
#' models have identical autocorrelations, partial autocorrelations and
#' Gaussian copula likelihoods.
#'
#' @param object an [arfactcop-class] model with specified coefficients, or
#'   an [arfactcopfit-class] object (its fitted model is converted; the
#'   standard errors are not carried over).
#'
#' @return An [armacop-class] object with `ar` coefficients
#'   \eqn{\prod_i (1 - \lambda_i B)} and `ma` coefficients from the
#'   factorisation (in the [stats::arima()] sign convention).
#' @export
#'
#' @examples
#' g <- arfactcop(lambda = c(0.98, 0.7), w = c(0.3, 0.4))
#' a <- arfact2arma(g)
#' a
#' # identical partial autocorrelations
#' max(abs(ktau(g, 20) - ktau(a, 20)))
arfact2arma <- function(object) {
  if (is(object, "arfactcopfit")) object <- object@model
  if (!is(object, "arfactcop")) stop("'object' must be an arfactcop model or fit")
  if (!gc_is_specified(object)) stop("the model has unspecified coefficients; fit it to data first")
  lam <- object@lambda
  w <- object@w
  m <- object@m
  w0 <- if (object@noise) max(1 - sum(w), 0) else 0
  # AR polynomial prod (1 - lambda_i B)
  arpoly <- 1
  for (l in lam) arpoly <- poly_mult(arpoly, c(1, -l))
  ar <- -arpoly[-1L]
  q <- if (w0 > 1e-12) m else m - 1L
  if (q == 0L) return(armacop(ar = ar))
  # autocovariances of the MA(q) process y = phi(B) x
  rho <- c(1, afc_acf(lam, w, m + q))
  gx <- function(h) rho[abs(h) + 1L]
  phi <- arpoly # phi_0 = 1, phi_i
  gam <- vapply(0:q, function(k) {
    s <- 0
    for (i in 0:m) for (j in 0:m) s <- s + phi[i + 1L] * phi[j + 1L] * gx(k + j - i)
    s
  }, numeric(1))
  # roots of P(z) = sum_{k=-q}^{q} gam_|k| z^{q+k}; keep the q roots outside the unit circle
  cf <- c(rev(gam[-1L]), gam[1L], gam[-1L]) # ascending powers z^0 ... z^{2q}
  rts <- polyroot(cf)
  out <- rts[Mod(rts) > 1]
  if (length(out) != q) {
    stop("spectral factorisation failed: the moving-average polynomial has roots on the unit circle")
  }
  theta <- 1
  for (r in out) theta <- poly_mult(theta, c(1, -1 / r))
  ma <- Re(theta[-1L])
  armacop(ar = ar, ma = ma)
}

setAs("arfactcop", "armacop", function(from) arfact2arma(from))
