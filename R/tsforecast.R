# Out-of-sample forecasting: one-step predictive distributions at every origin of a holdout, with optional
# periodic re-estimation (tsroll), and the class that holds them (tsforecast). Scoring is in tsscore.R.

#' Class of out-of-sample forecasts
#'
#' Output of [tsroll()]: for every forecast origin the one-step-ahead predictive distribution of the next
#' observation, summarised by its quantiles, the probability integral transform (PIT) and the log predictive
#' density at the observed value, and the value-at-risk and expected shortfall at the requested tail
#' probabilities.
#'
#' @slot index positions in the data of the observations forecast.
#' @slot observed the observations forecast (on the scale of the data).
#' @slot levels the probability levels of the columns of `quantile`.
#' @slot quantile matrix (one row per forecast, one column per level) of predictive quantiles.
#' @slot pit the PIT of each observation under its predictive distribution (iid uniform for a correct model).
#' @slot logdens the log predictive density at each observation (for a model with a margin, including the
#'   log density of the margin; for an empirical margin its kernel density estimate).
#' @slot p the tail probabilities of `var` and `es`.
#' @slot var matrix of lower-tail value-at-risk (the predictive quantile at `p`), one column per `p`.
#' @slot es matrix of lower-tail expected shortfall, one column per `p`.
#' @slot blocks a data frame with one row per block of forecasts: the first and last origin, the window of the
#'   data used to estimate the model and the status of the (re-)estimation (`"given"` for the model of the
#'   call, `"ok"`, or `"failed"` when the previous model was kept).
#' @slot estimates matrix of the coefficients of the model of each block (one row per block); empty if the
#'   model has no coefficient vector.
#' @slot window,refit,width the settings of the call.
#' @slot model the model of the last block.
#'
#' @exportClass tsforecast
setClass("tsforecast", slots = c(
  index = "integer", observed = "numeric", levels = "numeric", quantile = "matrix", pit = "numeric",
  logdens = "numeric", p = "numeric", var = "matrix", es = "matrix", blocks = "data.frame",
  estimates = "matrix", window = "character", refit = "numeric", width = "numeric", model = "ANY"
))

#' @rdname tsforecast-class
#' @param object a [tsforecast-class] object.
#' @export
setMethod("show", "tsforecast", function(object) {
  n <- length(object@index)
  cat("object class: tsforecast\n")
  cat(n, " forecasts (positions ", min(object@index), " to ", max(object@index), "), ", length(object@levels),
      " quantile levels\n", sep = "")
  rf <- if (is.finite(object@refit)) sprintf("re-estimated every %d observations on %s windows%s", as.integer(object@refit),
                                             object@window, if (object@window == "rolling") sprintf(" of width %d", as.integer(object@width)) else "")
        else "parameters held fixed"
  cat(nrow(object@blocks), " block(s): ", rf, "\n", sep = "")
  cat("mean log predictive density: ", format(mean(object@logdens), digits = 5), "\n", sep = "")
  if (length(object@p) > 0L) {
    br <- colSums(object@observed < object@var)
    cat("lower-tail VaR breaches: ", paste0(format(object@p), ": ", br, " (", format(object@p * n, digits = 3), " expected)", collapse = "; "), "\n", sep = "")
  }
  invisible(object)
})

# ---- the model behind an object ---------------------------------------------------------------------------

# Specified model, its kind and default data from a model or a fitted object.
tr_unpack <- function(object, data) {
  if (is(object, "tsmodfit")) { model <- object@model; if (missing(data) || is.null(data)) data <- tm_fit_data(object) }
  else if (is(object, "tsmod")) model <- object
  else if (is(object, "sdvinecopfit") || is(object, "vtscopfit") || is(object, "armacopfit") || is(object, "arfactcopfit")) {
    model <- object@model
    if (missing(data) || is.null(data)) data <- object@data
  } else if (is(object, "sdvinecop") || is(object, "vtscop") || is(object, "armacop") || is(object, "arfactcop")) model <- object
  else stop("'object' must be a tsmod, sdvinecop, vtscop, armacop or arfactcop model or a fitted object of one of these classes", call. = FALSE)
  if (missing(data) || is.null(data)) stop("'data' is needed", call. = FALSE)
  if (is(model, "sdvinecop") && model@mode != "tied") stop("free-mode d-vine models are not supported", call. = FALSE)
  list(model = model, data = as.numeric(data), is_tsmod = is(model, "tsmod"))
}

tr_copula <- function(model) if (is(model, "tsmod")) model@copula else model

# ---- the pass at fixed parameters --------------------------------------------------------------------------

# Predictive quantiles `lev` (u scale), PIT and log density of u[t] for each origin t, given the previous values.
tr_pass <- function(cop, u, origins, lev, cores = 1L) {
  if (is(cop, "sdvinecop")) {
    pairs <- dv_resolve(cop)
    if (is.null(pairs)) stop("the model is infeasible at its coefficients", call. = FALSE)
    e <- if (isTRUE(getOption("tscop.cpp", TRUE))) dv_encode(pairs, ast = TRUE, nl = TRUE) else NULL
    if (!is.null(e)) {
      chunk <- function(o) dv_holdout_cpp(u, o, lev, e$fam, e$rot, e$par, e$ivt, e$bfam, e$brot, e$bpar, e$d1, e$d2, e$nl, e$k1, e$k2)
      return(tr_chunks(origins, chunk, cores))
    }
  }
  # any copula process with a predictive distribution: one call of predict() per origin
  one <- function(t) {
    if (t == 1L) return(list(pit = u[1L], logdens = 0, quantile = lev))
    past <- u[seq_len(t - 1L)]
    list(pit = as.numeric(predict(cop, past, type = "df", x = u[t])),
         logdens = log(as.numeric(predict(cop, past, type = "dens", x = u[t]))),
         quantile = as.numeric(predict(cop, past, type = "qf", x = lev)))
  }
  res <- if (cores > 1L) parallel::mclapply(origins, one, mc.cores = cores) else lapply(origins, one)
  list(pit = vapply(res, `[[`, 0, "pit"), logdens = vapply(res, `[[`, 0, "logdens"),
       quantile = do.call(rbind, lapply(res, `[[`, "quantile")))
}

# Run a compiled holdout pass over chunks of the origins (in parallel with cores > 1) and bind the results.
tr_chunks <- function(origins, f, cores) {
  if (cores > 1L && length(origins) >= 2L * cores && .Platform$OS.type != "windows") {
    parts <- split(origins, cut(seq_along(origins), cores, labels = FALSE))
    res <- parallel::mclapply(parts, function(o) f(as.integer(o)), mc.cores = cores)
    if (any(vapply(res, inherits, TRUE, "try-error"))) stop("a parallel holdout pass failed", call. = FALSE)
    list(pit = unlist(lapply(res, `[[`, "pit")), logdens = unlist(lapply(res, `[[`, "logdens")),
         quantile = do.call(rbind, lapply(res, `[[`, "quantile")))
  } else {
    f(as.integer(origins))
  }
}

# Expected-shortfall quadrature nodes for tail probability p: ES = sum_k w_k q(s_k) exp(-v_k) with
# s_k = p exp(-v_k) (as tsrisk()).
tr_es_nodes <- function(p, nodes) {
  g <- gl_nodes(nodes)
  vmax <- min(40, log(p / 1e-14))
  list(s = p * exp(-g$s * vmax), wts = g$w * vmax * exp(-g$s * vmax))
}

# One block: forecasts of the observations at positions `origins` with the model `model`, given all earlier data.
tr_block <- function(model, x, origins, levels, p, es_nodes, cores) {
  is_t <- is(model, "tsmod")
  cop <- tr_copula(model)
  mg <- if (is_t) model@margin else NULL
  u <- if (is_t) tm_pits(mg, x[seq_len(max(origins))]) else x[seq_len(max(origins))]
  en <- lapply(p, tr_es_nodes, nodes = es_nodes)
  lev <- c(levels, p, unlist(lapply(en, `[[`, "s")))
  r <- tr_pass(cop, u, origins, lev, cores)
  Q <- r$quantile
  if (is_t) Q <- matrix(mqf(mg, as.numeric(Q)), nrow(Q), ncol(Q))
  nl <- length(levels); np <- length(p)
  out <- list(pit = r$pit, logdens = r$logdens + if (is_t) log(mdens(mg, x[origins])) else 0,
              quantile = Q[, seq_len(nl), drop = FALSE], var = Q[, nl + seq_len(np), drop = FALSE])
  out$es <- vapply(seq_len(np), function(i) {
    cols <- nl + np + (i - 1L) * es_nodes + seq_len(es_nodes)
    as.numeric(Q[, cols, drop = FALSE] %*% en[[i]]$wts)
  }, numeric(nrow(Q)))
  out$es <- matrix(out$es, nrow(Q), np)
  out
}

# Re-estimate the model on the data of a window; the previous coefficients are the starting values. NULL on failure.
tr_refit <- function(model, xw, fit_args) {
  tryCatch(suppressWarnings({
    if (is(model, "tsmod")) {
      init <- if (is(model@copula, "sdvinecop") || is(model@copula, "vtscop")) list(init = model@copula) else list()
      do.call(tsfit, c(list(xw, model), init, fit_args))@model
    } else if (is(model, "sdvinecop") || is(model, "vtscop")) {
      do.call(tsfit, c(list(xw, model, pseudo = FALSE, init = model), fit_args))@model
    } else {
      do.call(tsfit, c(list(xw, model, pseudo = FALSE), fit_args))@model
    }
  }), error = function(e) NULL)
}

#' Out-of-sample forecasts with periodic re-estimation
#'
#' For every forecast origin in a holdout, the one-step-ahead predictive distribution of the next observation
#' given all the earlier data, summarised for forecast evaluation (see [tsscore()]): predictive quantiles on a
#' grid of probability levels, the probability integral transform (PIT) and the log predictive density at the
#' observed value, and the value-at-risk and expected shortfall of the lower tail. With `refit = Inf` the
#' parameters stay at those of `object` (a fixed-parameter holdout, computed in one compiled pass for
#' [sdvinecop-class] copulas); with `refit = H` the model is re-estimated every `H` observations on an
#' expanding or rolling window ending just before the block, starting from the previous coefficients, while the
#' forecasts always condition on all the observed data before each origin.
#'
#' The model of the first block is the model of `object`; later blocks use the model re-estimated on the window
#' before the block (the family and structure of `object`, with its coefficients as starting values). A
#' [tsmod-class] model has its margin re-estimated as well (an empirical margin is recomputed from the window);
#' a copula model is fitted to its data on the copula scale (`pseudo = FALSE`). When a re-estimation fails the
#' previous model is kept and the block is marked `"failed"`.
#'
#' The expected shortfall of the lower tail is \eqn{p^{-1}\int_0^p F^{-1}(s)\,ds}, computed as in [tsrisk()] by
#' Gauss-Legendre quadrature of the predictive quantile function after the substitution \eqn{s = p e^{-v}}
#' (`es_nodes` nodes). For losses, forecast the negated series.
#'
#' @param object a [tsmod-class], [sdvinecop-class], [vtscop-class], [armacop-class] or [arfactcop-class] model
#'   with known coefficients, or a fitted object of one of these classes.
#' @param data the complete series: observations on their own scale for a [tsmod-class] model, copula data in
#'   (0, 1) for a copula model. For a fitted object it defaults to the data of the fit.
#' @param start position in `data` of the first observation to forecast (at least 2); the forecasts are for
#'   `start, ..., length(data)`.
#' @param refit re-estimation interval in observations, or `Inf` (default) for fixed parameters.
#' @param window `"expanding"` (default): all the data before the block; `"rolling"`: the last `width`
#'   observations before the block.
#' @param width the width of a rolling window (default 1500, or the number of observations before `start` if
#'   smaller).
#' @param levels probability levels of the predictive quantiles (default `(1:99)/100`).
#' @param p tail probabilities of the value-at-risk and expected shortfall (default `c(0.01, 0.05)`). With `p = NULL`
#'   no risk measures are computed: the forecasts then have no value-at-risk or expected shortfall, [tsscore()] scores only
#'   the forecasts, and the computation is faster, because the expected shortfall needs many extra predictive quantiles.
#' @param es_nodes number of quadrature nodes of the expected shortfall (default 32).
#' @param fit_args further arguments passed to [tsfit()] in the re-estimations (for example
#'   `list(method = "joint")` for a [tsmod-class] model).
#' @param cores number of processes (not on Windows). The compiled pass is split over the origins; the
#'   re-estimations are sequential because each starts from the previous one.
#'
#' @return A [tsforecast-class] object.
#' @seealso [tsscore()] for the scores and tests, [tsrisk()] for the risk measures at the end of the data.
#' @export
#' @examples
#' \donttest{
#' set.seed(1)
#' m <- tsmod(marg("st", c(mu = 0, sigma = 0.01, df = 5)),
#'            sdvinecop(arfactcop(lambda = 0.95, w = 0.3), tail = pc("ivt", base = "joe"),
#'                      vt = c(0.6, 0.45), maxlag = 12, tautol = 1e-3))
#' x <- tssim(m, 1000)
#' fc <- tsroll(m, x, start = 801)
#' fc
#' tsscore(fc)
#' }
tsroll <- function(object, data, start, refit = Inf, window = c("expanding", "rolling"), width = 1500L,
                   levels = (1:99) / 100, p = c(0.01, 0.05), es_nodes = 32L, fit_args = list(), cores = 1L) {
  window <- match.arg(window)
  up <- tr_unpack(object, if (missing(data)) NULL else data)
  model <- up$model
  x <- up$data
  n <- length(x)
  if (anyNA(x) || any(!is.finite(x))) stop("'data' must be complete and finite", call. = FALSE)
  if (!up$is_tsmod && any(x <= 0 | x >= 1)) stop("copula data must lie in (0, 1)", call. = FALSE)
  if (missing(start) || !is.numeric(start) || length(start) != 1L || start < 2 || start > n || start != round(start)) {
    stop("'start' must be a whole number between 2 and the length of the data", call. = FALSE)
  }
  if (!is.numeric(levels) || anyNA(levels) || any(levels <= 0 | levels >= 1)) stop("'levels' must be probabilities strictly inside (0, 1)", call. = FALSE)
  if (is.null(p)) p <- numeric(0)
  if (!is.numeric(p) || anyNA(p) || any(p <= 0 | p >= 1)) stop("'p' must be probabilities strictly inside (0, 1), or NULL", call. = FALSE)
  if (!is.numeric(refit) || length(refit) != 1L || refit < 1 || (is.finite(refit) && refit != round(refit))) {
    stop("'refit' must be a positive whole number or Inf", call. = FALSE)
  }
  es_nodes <- as.integer(es_nodes)
  width <- as.integer(min(width, start - 1L))
  origins <- start:n
  blocks <- if (is.finite(refit)) split(origins, ceiling(seq_along(origins) / refit)) else list(origins)
  res <- vector("list", length(blocks))
  meta <- data.frame(block = seq_along(blocks), from = NA_integer_, to = NA_integer_, window_from = NA_integer_,
                     window_to = NA_integer_, status = "given", stringsAsFactors = FALSE)
  ests <- list()
  cur <- model
  for (b in seq_along(blocks)) {
    o <- blocks[[b]]
    t0 <- o[1L]
    a <- if (window == "expanding") 1L else max(1L, t0 - width)
    meta$from[b] <- t0; meta$to[b] <- o[length(o)]; meta$window_from[b] <- a; meta$window_to[b] <- t0 - 1L
    if (b > 1L) {
      nm <- tr_refit(cur, x[a:(t0 - 1L)], fit_args)
      if (is.null(nm)) meta$status[b] <- "failed" else { cur <- nm; meta$status[b] <- "ok" }
    }
    res[[b]] <- tr_block(cur, x, o, levels, p, es_nodes, cores)
    ests[[b]] <- tryCatch(coef(cur), error = function(e) numeric(0))
  }
  bind <- function(f) do.call(rbind, lapply(res, `[[`, f))
  est <- if (all(vapply(ests, length, 1L) == length(ests[[1L]])) && length(ests[[1L]]) > 0L) do.call(rbind, ests) else matrix(numeric(0), length(ests), 0L)
  new("tsforecast", index = as.integer(origins), observed = x[origins], levels = levels,
      quantile = bind("quantile"), pit = unlist(lapply(res, `[[`, "pit")), logdens = unlist(lapply(res, `[[`, "logdens")),
      p = p, var = bind("var"), es = bind("es"), blocks = meta, estimates = est, window = window,
      refit = as.numeric(refit), width = as.numeric(width), model = cur)
}
