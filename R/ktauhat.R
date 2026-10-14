# Empirical Kendall partial autocorrelations of data given a fitted or specified model.

#' Empirical Kendall partial autocorrelations of data given a model
#'
#' The "generalised lag" diagnostic: for lag 1 the pairs \eqn{(u_t, u_{t+1})}; at
#' each further lag the pairs are conditioned on the observations in between
#' through the h-functions of the model's pair copulas, so that the pairs at lag
#' \eqn{k} are the arguments the model feeds to its lag-\eqn{k} pair copula. The
#' Kendall's tau of these pairs estimates the Kendall partial autocorrelation at
#' lag \eqn{k} if the pair copulas at lags below \eqn{k} are right, whatever the
#' family used at lag \eqn{k}. Comparing it with the model's own
#' [ktau()] shows where the fitted dependence structure departs from the data.
#'
#' For models with inverse-v-transformed pair copulas the pairs are v-transformed
#' (the dependence of the u-scale pairs is of cross shape and has tau zero): the
#' tau is that of the base copula pairs, the scale of the kpacf, as for
#' [ktau()]. For a [vtscop-class] model the data are v-transformed and the
#' diagnostic is that of the base process; for a [tsmod-class] model the data
#' (on the original scale) are first transformed by the margin.
#'
#' With `B > 0` the diagnostic is repeated for `B` series simulated from the model
#' (the parameters are taken as given, so parameter uncertainty is ignored) and the
#' pointwise `level` interval of the simulated values is returned in `lo` and `hi`: if the
#' data's values fall inside, they are what the model typically produces.
#'
#' @param object a [armacop-class], [arfactcop-class], [sdvinecop-class],
#'   [vtscop-class] or [tsmod-class] model with known coefficients, or a fitted
#'   model object.
#' @param data the data: copula data in (0, 1) without missing values, or for a
#'   [tsmod-class] model observations on the original scale. For a fitted object
#'   the data used in the fit by default.
#' @param lagmax the largest lag (default 20; at most the number of observations
#'   minus 12).
#' @param pairs the number of lags (at most 9, default 0) for which the pairs
#'   themselves are returned, as an attribute `"pairs"`: a list of two-column
#'   matrices on the copula scale (the first lags of the legacy generalised lag
#'   plots; for a [vtscop-class] model the pairs of the v-transformed data).
#' @param B the number of simulated series for the interval (default 0: none).
#' @param level the coverage of the interval (default 0.95).
#' @param cores the number of processes for the simulated series (default 1; not on Windows). The
#'   result does not depend on it: each replicate has its own seed, drawn from the random numbers of the caller.
#' @param ... further arguments passed to methods.
#'
#' @return A data frame with the columns `lag`, `tau_data` (the empirical
#'   Kendall partial autocorrelation), `tau_model` (the model's kpacf in Kendall's
#'   tau, as [ktau()]), `n` (the number of pairs) and, with `B > 0`, `lo` and `hi`.
#' @export
#' @examples
#' set.seed(1)
#' m <- sdvinecop(arfactcop(lambda = 0.9, w = 0.5), pairs = list(pc("gumbel")))
#' u <- tssim(m, 1000)
#' ktauhat(m, u, lagmax = 8)
#' # the model fitted to other data: the empirical values depart from its kpacf
#' ktauhat(sdvinecop(arfactcop(lambda = 0.5, w = 0.5), pairs = list(pc("gumbel"))), u, lagmax = 8, B = 50)
setGeneric("ktauhat", function(object, data, lagmax = 20, ...) standardGeneric("ktauhat"),
           signature = "object")

# Kendall's tau-b, O(m log m) in C++ (stats::cor is O(m^2)) unless options(tscop.cpp = FALSE)
kendall_tau <- function(x, y) {
  if (isTRUE(getOption("tscop.cpp", TRUE))) kendall_tau_cpp(as.numeric(x), as.numeric(y))
  else stats::cor(x, y, method = "kendall")
}

# R reference of the level recursion on the u scale: Kendall's tau of the level pairs (of the
# v-transformed pairs for IVT pair copulas) and the pairs themselves for the first `keep` levels
dv_level_tau_r <- function(pairs, u, keep = 0L) {
  n <- length(u)
  K <- min(length(pairs), n - 2L)
  A <- u[-n]
  B <- u[-1L]
  tau <- rep(NA_real_, length(pairs))
  kept <- list()
  for (j in seq_len(K)) {
    pr <- pairs[[j]]
    tau[j] <- if (dv_is_ivt(pr)) kendall_tau(dv_vlin(A, pr$delta[1L]), dv_vlin(B, pr$delta[2L])) else kendall_tau(A, B)
    if (j <= keep) kept[[j]] <- cbind(A, B)
    if (j < K) {
      m <- length(A)
      Anew <- dv_h(pr, A[-m], B[-m], 2L)
      Bnew <- dv_h(pr, A[-1L], B[-1L], 1L)
      A <- Anew
      B <- Bnew
    }
  }
  list(tau = tau, pairs = kept)
}

# Level taus of the pair copulas `pairs` at the data u: compiled where every pair copula is supported
dv_level_tau <- function(pairs, u) {
  if (isTRUE(getOption("tscop.cpp", TRUE)) && length(pairs) > 0L) {
    e <- dv_encode(pairs, ast = TRUE)
    if (!is.null(e)) return(dv_level_tau_cpp(as.numeric(u), e$fam, e$rot, e$par, e$ivt, e$bfam, e$brot, e$bpar, e$d1, e$d2))
  }
  dv_level_tau_r(pairs, u)$tau
}

# Gaussian pair copulas with the partial autocorrelations of a Gaussian process
ktauhat_pairs_gauss <- function(G, lagmax) {
  lapply(gc_pacf(G, lagmax), function(a) list(family = "gauss", rotation = 0, parameters = a))
}

# Pair copulas of a d-vine model at lags 1, ..., lagmax (beyond the effective lag K with the tail family)
ktauhat_pairs_dv <- function(model, lagmax) {
  m <- model
  m@K <- as.integer(max(lagmax, length(model@pairs), 1L))
  pr <- dv_resolve(m)
  if (is.null(pr)) stop("the model is infeasible at its current parameter values")
  pr[seq_len(lagmax)]
}

ktauhat_check <- function(u, lagmax, pairs) {
  u <- as.numeric(u)
  if (anyNA(u) || any(!is.finite(u)) || any(u <= 0 | u >= 1)) stop("'data' must be copula data in (0, 1) without missing values", call. = FALSE)
  if (!is.numeric(lagmax) || length(lagmax) != 1L || lagmax < 1 || lagmax != round(lagmax)) stop("'lagmax' must be a positive whole number", call. = FALSE)
  L <- as.integer(min(lagmax, length(u) - 12L))
  if (L < 1L) stop("too few observations", call. = FALSE)
  keep <- if (is.null(pairs)) 0L else as.integer(pairs)
  if (length(keep) != 1L || is.na(keep) || keep < 0L) stop("'pairs' must be a non-negative whole number", call. = FALSE)
  if (keep > 9L) {
    warning("'pairs' is capped at 9")
    keep <- 9L
  }
  list(u = u, L = L, keep = min(keep, L))
}

# Common worker: `pairs` the pair copulas, `u` the (possibly v-transformed) data on the scale of the pair
# copulas, `tau_model` the model's kpacf, `sim1()` one simulated series on the scale of `u`
ktauhat_core <- function(pairs, u, tau_model, L, keep, B, level, sim1, cores) {
  tau <- dv_level_tau(pairs, u)
  out <- data.frame(lag = seq_len(L), tau_data = tau[seq_len(L)], tau_model = tau_model[seq_len(L)],
                    n = length(u) - seq_len(L))
  if (B > 0L) {
    if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1) stop("'level' must be in (0, 1)", call. = FALSE)
    # one seed per replicate, drawn from the caller's random numbers: the result does not depend on `cores`
    seeds <- sample.int(.Machine$integer.max, B)
    one <- function(b) {
      had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
      old <- if (had) get(".Random.seed", envir = globalenv()) else NULL
      on.exit(if (had) assign(".Random.seed", old, envir = globalenv()))
      set.seed(seeds[b])
      dv_level_tau(pairs, sim1())[seq_len(L)]
    }
    res <- if (cores > 1L && .Platform$OS.type != "windows") parallel::mclapply(seq_len(B), one, mc.cores = cores) else lapply(seq_len(B), one)
    ts <- matrix(unlist(res), nrow = L)
    out$lo <- apply(ts, 1L, stats::quantile, probs = (1 - level) / 2, na.rm = TRUE, names = FALSE)
    out$hi <- apply(ts, 1L, stats::quantile, probs = 1 - (1 - level) / 2, na.rm = TRUE, names = FALSE)
    attr(out, "level") <- level
  }
  if (keep > 0L) attr(out, "pairs") <- dv_level_tau_r(pairs[seq_len(keep)], u, keep)$pairs
  out
}

ktauhat_known <- function(model) {
  if (is(model, "gausscop") && !gc_is_specified(model)) stop("the model has unspecified coefficients; fit it to data first", call. = FALSE)
}

#' @describeIn ktauhat Gaussian ARMA copula: Gaussian pair copulas with the partial autocorrelations of the process.
#' @export
setMethod("ktauhat", "armacop", function(object, data, lagmax = 20, pairs = 0L, B = 0L, level = 0.95, cores = 1L, ...) {
  if (missing(data)) stop("'data' is required", call. = FALSE)
  ktauhat_known(object)
  ck <- ktauhat_check(data, lagmax, pairs)
  n <- length(ck$u)
  ktauhat_core(ktauhat_pairs_gauss(object, ck$L), ck$u, ktau(object, ck$L), ck$L, ck$keep, B, level,
               function() tssim(object, n), cores)
})

#' @describeIn ktauhat AR-factor copula (as for the ARMA copula).
#' @export
setMethod("ktauhat", "arfactcop", function(object, data, lagmax = 20, pairs = 0L, B = 0L, level = 0.95, cores = 1L, ...) {
  if (missing(data)) stop("'data' is required", call. = FALSE)
  ktauhat_known(object)
  ck <- ktauhat_check(data, lagmax, pairs)
  n <- length(ck$u)
  ktauhat_core(ktauhat_pairs_gauss(object, ck$L), ck$u, ktau(object, ck$L), ck$L, ck$keep, B, level,
               function() tssim(object, n), cores)
})

#' @describeIn ktauhat Stationary d-vine copula: the model's pair copulas at lags `1, ..., lagmax` (the tail
#'   family beyond the effective lag).
#' @export
setMethod("ktauhat", "sdvinecop", function(object, data, lagmax = 20, pairs = 0L, B = 0L, level = 0.95, cores = 1L, ...) {
  if (missing(data)) stop("'data' is required", call. = FALSE)
  ck <- ktauhat_check(data, lagmax, pairs)
  n <- length(ck$u)
  pc_list <- ktauhat_pairs_dv(object, ck$L)
  tm <- if (object@mode == "tied") dv_kpacf_tau(object@kpacf, ck$L) else ktau(object, ck$L)
  sim_pairs <- dv_pairs_or_stop(object)
  ktauhat_core(pc_list, ck$u, tm, ck$L, ck$keep, B, level,
               function() as.numeric(dv_simulate(sim_pairs, matrix(stats::runif(n), n, 1L))), cores)
})

#' @describeIn ktauhat Process-level v-transform copula: the diagnostic of the base process for the
#'   v-transformed data.
#' @export
setMethod("ktauhat", "vtscop", function(object, data, lagmax = 20, pairs = 0L, B = 0L, level = 0.95, cores = 1L, ...) {
  if (missing(data)) stop("'data' is required", call. = FALSE)
  if (!vt_base_specified(object)) stop("the model has unspecified coefficients; fit it to data first", call. = FALSE)
  ck <- ktauhat_check(data, lagmax, pairs)
  n <- length(ck$u)
  v <- vt_transform_data(object@vt, ck$u)
  b <- object@base
  pc_list <- if (is(b, "arfactcop")) ktauhat_pairs_gauss(b, ck$L) else ktauhat_pairs_dv(b, ck$L)
  tm <- if (is(b, "arfactcop")) ktau(b, ck$L) else dv_kpacf_tau(b@kpacf, ck$L)
  ktauhat_core(pc_list, v, tm, ck$L, ck$keep, B, level,
               function() vt_transform_data(object@vt, vt_simulate(object, n, 1L)[, 1L]), cores)
})

#' @describeIn ktauhat Full time series model: the data (original scale) are transformed by the margin and the
#'   diagnostic is that of the copula process.
#' @export
setMethod("ktauhat", "tsmod", function(object, data, lagmax = 20, pairs = 0L, B = 0L, level = 0.95, cores = 1L, ...) {
  if (missing(data)) stop("'data' is required", call. = FALSE)
  if (tm_auto(object@margin)) stop("the margin is chosen automatically when the model is fitted; fit it to data first", call. = FALSE)
  tm_check_data(data)
  ktauhat(object@copula, tm_pits(object@margin, data), lagmax, pairs = pairs, B = B, level = level, cores = cores)
})

#' @describeIn ktauhat Fitted model objects: the model and, by default, the data of the fit.
#' @export
setMethod("ktauhat", "armacopfit", function(object, data, lagmax = 20, ...) {
  if (missing(data)) data <- object@data
  ktauhat(object@model, data, lagmax, ...)
})
#' @rdname ktauhat
#' @export
setMethod("ktauhat", "arfactcopfit", function(object, data, lagmax = 20, ...) {
  if (missing(data)) data <- object@data
  ktauhat(object@model, data, lagmax, ...)
})
#' @rdname ktauhat
#' @export
setMethod("ktauhat", "sdvinecopfit", function(object, data, lagmax = 20, ...) {
  if (missing(data)) data <- object@data
  ktauhat(object@model, data, lagmax, ...)
})
#' @rdname ktauhat
#' @export
setMethod("ktauhat", "vtscopfit", function(object, data, lagmax = 20, ...) {
  if (missing(data)) data <- object@data
  ktauhat(object@model, data, lagmax, ...)
})
#' @rdname ktauhat
#' @export
setMethod("ktauhat", "tsmodfit", function(object, data, lagmax = 20, ...) {
  if (missing(data)) data <- object@data
  ktauhat(object@model, data, lagmax, ...)
})
