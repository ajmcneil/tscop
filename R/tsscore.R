# Scoring and backtesting of out-of-sample forecasts (see tsroll() in tsforecast.R).

# Newey-West (Bartlett) long-run variance of a series, mean removed.
ts_lrv <- function(d, lag = NULL) {
  n <- length(d)
  if (is.null(lag)) lag <- floor(n^(1 / 3))
  d <- d - mean(d)
  g0 <- sum(d * d) / n
  if (lag < 1L) return(g0)
  g <- vapply(seq_len(lag), function(k) sum(d[-seq_len(k)] * d[seq_len(n - k)]) / n, 0)
  g0 + 2 * sum((1 - seq_len(lag) / (lag + 1)) * g)
}

# z-test of mean zero for a loss differential (or other martingale-type series)
ts_mean_test <- function(d, lag = NULL) {
  n <- length(d)
  v <- ts_lrv(d, lag)
  z <- if (is.finite(v) && v > 0) mean(d) / sqrt(v / n) else NA_real_   # NA for non-finite differences (for example an infinite log score)
  c(mean = mean(d), stat = z, p.value = 2 * stats::pnorm(-abs(z)))
}

# Kupiec and Christoffersen likelihood-ratio tests for breach indicators I with nominal rate p
ts_var_tests <- function(I, p) {
  n <- length(I); x <- sum(I)
  ll <- function(pi_, x, n) if (pi_ <= 0 || pi_ >= 1) (if ((pi_ <= 0 && x == 0) || (pi_ >= 1 && x == n)) 0 else -Inf) else x * log(pi_) + (n - x) * log1p(-pi_)
  lr_uc <- -2 * (ll(p, x, n) - ll(x / n, x, n))
  a <- I[-n]; b <- I[-1L]
  n00 <- sum(a == 0 & b == 0); n01 <- sum(a == 0 & b == 1); n10 <- sum(a == 1 & b == 0); n11 <- sum(a == 1 & b == 1)
  pi01 <- if (n00 + n01 > 0) n01 / (n00 + n01) else 0
  pi11 <- if (n10 + n11 > 0) n11 / (n10 + n11) else 0
  pi_all <- (n01 + n11) / max(1, n00 + n01 + n10 + n11)
  l1 <- ll(pi01, n01, n00 + n01) + ll(pi11, n11, n10 + n11)
  l0 <- ll(pi_all, n01 + n11, n00 + n01 + n10 + n11)
  lr_ind <- max(0, -2 * (l0 - l1))
  c(breaches = x, expected = n * p, rate = x / n, kupiec = lr_uc, kupiec.p = stats::pchisq(lr_uc, 1, lower.tail = FALSE),
    independence = lr_ind, independence.p = stats::pchisq(lr_ind, 1, lower.tail = FALSE),
    conditional = lr_uc + lr_ind, conditional.p = stats::pchisq(lr_uc + lr_ind, 2, lower.tail = FALSE))
}

# Per-time scores of one forecast: pinball losses by level, AQS, FZ0, etc.
ts_unpack_fc <- function(f, observed, levels, p) {
  if (is(f, "tsforecast")) {
    # value-at-risk and expected shortfall at the requested probabilities p (a subset of those of the forecast)
    var <- f@var; es <- f@es; pp <- f@p
    if (length(p) > 0L) {
      j <- vapply(p, function(a) { k <- which(abs(f@p - a) < 1e-9); if (length(k) == 1L) k else NA_integer_ }, 1L)
      if (anyNA(j)) stop("the forecasts contain the value-at-risk and expected shortfall only at p = ", paste(f@p, collapse = ", "), call. = FALSE)
      var <- var[, j, drop = FALSE]; es <- es[, j, drop = FALSE]; pp <- p
    }
    list(y = f@observed, Q = f@quantile, lev = f@levels, logdens = f@logdens, p = pp, var = var, es = es)
  } else {
    Q <- as.matrix(f)
    if (is.null(observed)) stop("'observed' is needed for forecasts given as quantile matrices", call. = FALSE)
    if (is.null(levels)) stop("'levels' is needed for forecasts given as quantile matrices", call. = FALSE)
    if (ncol(Q) != length(levels)) stop("the number of columns of the quantile matrix must equal the number of levels", call. = FALSE)
    if (nrow(Q) != length(observed)) stop("the quantile matrix must have one row per observation", call. = FALSE)
    V <- matrix(NA_real_, nrow(Q), length(p))
    for (i in seq_along(p)) {
      j <- which(abs(levels - p[i]) < 1e-9)
      if (length(j) == 1L) V[, i] <- Q[, j]
    }
    list(y = as.numeric(observed), Q = Q, lev = as.numeric(levels), logdens = rep(NA_real_, nrow(Q)), p = p, var = V,
         es = matrix(NA_real_, nrow(Q), length(p)))
  }
}

# pinball losses (n x L)
ts_pinball <- function(y, Q, lev) {
  ind <- outer(y, rep(1, length(lev))) < Q
  (ind - rep(lev, each = length(y))) * (Q - y)
}

#' Scores and backtests of out-of-sample forecasts
#'
#' Evaluates one or several sets of one-step-ahead forecasts of the same observations: the average quantile
#' score (AQS) over a grid of probability levels together with its profile by level, the negative log score,
#' backtests of the value-at-risk (Kupiec and Christoffersen likelihood-ratio tests of the breach indicators),
#' backtests of the expected shortfall (the Fissler-Ziegel loss FZ0 and the Acerbi-Szekely statistics Z1
#' and Z2), and Diebold-Mariano tests of equal predictive accuracy against the first forecast, with
#' Newey-West standard errors.
#'
#' @details
#' The quantile score at level \eqn{\tau} is the pinball loss \eqn{(1\{y<q\}-\tau)(q-y)}, averaged over the
#' forecasts; AQS is its weighted average over the levels (a lower AQS is better) and measures the whole
#' predictive distribution, weighting the levels as `weights` says (`"uniform"`, or non-negative numbers, or a
#' function of the level, normalised to sum to one). The VaR tests use the breach indicators
#' \eqn{1\{y<\mathrm{VaR}_p\}}; `kupiec` tests the breach rate, `independence` the absence of clustering of
#' breaches and `conditional` both. The ES statistics are for the lower tail with negative expected shortfall:
#' `fz0` is the mean Fissler-Ziegel loss of the pair (VaR, ES) (lower is better); `z1` is the mean of
#' \eqn{y/\mathrm{ES}-1} over breaches and `z2` the mean of
#' \eqn{y1\{y<\mathrm{VaR}\}/(p\,\mathrm{ES})-1}, both zero for correct forecasts and positive when the
#' expected shortfall is underestimated, with z-tests (Newey-West for `z2`). The Diebold-Mariano
#' tests compare per-forecast losses (AQS, the negative log score and FZ0 at each `p`): the statistic is
#' the mean of the loss of the first forecast minus the loss of the other, so positive values mean the other
#' forecast is better.
#'
#' All scores in the `summary` are losses: smaller (or more negative) is better. The log score is therefore
#' reported as the negative log score \eqn{-\log f(y)}, averaged over the forecasts and named `neglogscore`.
#'
#' @param x a [tsforecast-class] object, or a matrix of predictive quantiles with one row per observation and
#'   one column per level.
#' @param ... further forecasts of the same observations (each a [tsforecast-class] object or quantile matrix),
#'   named for the output (an unnamed argument that is a plain object name is labelled with it).
#' @param observed the observations, for forecasts given as matrices.
#' @param levels the probability levels of the columns, for forecasts given as matrices. The scoring uses the
#'   levels of the first forecast, which all must share.
#' @param p tail probabilities of the VaR and ES backtests (default the `p` of the first [tsforecast-class]
#'   object). For a matrix, VaR backtests need a column at level `p` and ES backtests are not available.
#' @param weights the level weights of the AQS (see Details).
#' @param lag truncation lag of the Newey-West variance (default \eqn{\lfloor n^{1/3}\rfloor}).
#'
#' @return An object of class `tsscore`, a list with `summary` (a data frame with a row per forecast: AQS, negative
#'   log score and FZ0 by `p`), `profile` (matrix of the mean quantile score by level and forecast), `var` and
#'   `es` (lists by `p` of data frames of the tests per forecast), `dm` (data frame of Diebold-Mariano
#'   tests against the first forecast), `levels`, `weights` and `n`. It has `print` and `plot` methods; the
#'   plot shows the profile of quantile scores by level, or with `relative = TRUE` the difference from
#'   the first forecast.
#' @seealso [tsroll()]
#' @export
#' @examples
#' \donttest{
#' set.seed(1)
#' m <- tsmod(marg("st", c(mu = 0, sigma = 0.01, df = 5)),
#'            sdvinecop(arfactcop(lambda = 0.95, w = 0.3), tail = pc("ivt", base = "joe"),
#'                      vt = c(0.6, 0.45), maxlag = 12, tautol = 1e-3))
#' x <- tssim(m, 1000)
#' tsscore(tsroll(m, x, start = 801))
#' }
tsscore <- function(x, ..., observed = NULL, levels = NULL, p = NULL, weights = "uniform", lag = NULL) {
  fcs <- c(list(x), list(...))
  # names: those given, else the argument when it is a plain name, else forecast1, forecast2, ...
  exprs <- c(list(substitute(x)), as.list(substitute(list(...)))[-1L])
  nm <- names(fcs)
  if (is.null(nm)) nm <- rep("", length(fcs))
  if (is.null(names(fcs)) || nm[1L] == "") nm[1L] <- if (is.name(exprs[[1L]])) as.character(exprs[[1L]]) else ""
  for (k in which(nm == "")) nm[k] <- if (is.name(exprs[[k]])) as.character(exprs[[k]]) else paste0("forecast", k)
  if (is.null(p)) p <- if (is(x, "tsforecast")) x@p else numeric(0)
  u <- lapply(fcs, ts_unpack_fc, observed = observed, levels = levels, p = p)
  y <- u[[1L]]$y; lev <- u[[1L]]$lev; n <- length(y)
  for (k in seq_along(u)) {
    if (length(u[[k]]$y) != n || any(abs(u[[k]]$y - y) > 1e-12 * (1 + abs(y)))) stop("all forecasts must be of the same observations", call. = FALSE)
    if (length(u[[k]]$lev) != length(lev) || any(abs(u[[k]]$lev - lev) > 1e-12)) stop("all forecasts must have the same levels", call. = FALSE)
  }
  w <- if (is.function(weights)) weights(lev) else if (is.character(weights) && identical(weights, "uniform")) rep(1, length(lev)) else weights
  if (!is.numeric(w) || length(w) != length(lev) || any(w < 0) || sum(w) <= 0) stop("'weights' must be \"uniform\", a function of the level, or non-negative numbers, one per level", call. = FALSE)
  w <- w / sum(w)
  K <- length(u); np <- length(p)
  pin <- lapply(u, function(f) ts_pinball(y, f$Q, lev))
  qs_t <- lapply(pin, function(P) as.numeric(P %*% w))
  profile <- vapply(pin, colMeans, numeric(length(lev)))
  dim(profile) <- c(length(lev), K); dimnames(profile) <- list(format(lev), nm)
  # FZ0 losses per time, forecast and p
  fz0_t <- lapply(u, function(f) {
    vapply(seq_len(np), function(i) {
      v <- f$var[, i]; e <- f$es[, i]
      if (anyNA(v) || anyNA(e) || any(e >= 0)) return(rep(NA_real_, n))
      (y <= v) * (y - v) / (p[i] * e) + v / e + log(-e) - 1
    }, numeric(n))
  })
  logs <- vapply(u, function(f) mean(f$logdens), 0)
  summ <- data.frame(forecast = nm, aqs = vapply(qs_t, mean, 0), neglogscore = -logs, stringsAsFactors = FALSE)
  for (i in seq_len(np)) summ[[paste0("fz0_", format(p[i]))]] <- vapply(fz0_t, function(m) mean(m[, i]), 0)
  var_tab <- es_tab <- list()
  for (i in seq_len(np)) {
    key <- format(p[i])
    var_tab[[key]] <- do.call(rbind, lapply(seq_len(K), function(k) {
      v <- u[[k]]$var[, i]
      if (anyNA(v)) return(data.frame(forecast = nm[k], t(stats::setNames(rep(NA_real_, 9), c("breaches", "expected", "rate", "kupiec", "kupiec.p", "independence", "independence.p", "conditional", "conditional.p")))))
      data.frame(forecast = nm[k], t(ts_var_tests(as.integer(y < v), p[i])))
    }))
    es_tab[[key]] <- do.call(rbind, lapply(seq_len(K), function(k) {
      v <- u[[k]]$var[, i]; e <- u[[k]]$es[, i]
      bad <- anyNA(v) || anyNA(e) || any(e >= 0)
      if (bad) return(data.frame(forecast = nm[k], fz0 = NA_real_, z1 = NA_real_, z1.p = NA_real_, z2 = NA_real_, z2.p = NA_real_))
      br <- y < v
      z1 <- if (any(br)) { d <- y[br] / e[br] - 1; c(mean(d), if (length(d) > 1L) mean(d) / (stats::sd(d) / sqrt(length(d))) else NA_real_) } else c(NA_real_, NA_real_)
      a <- ifelse(br, y / (p[i] * e), 0) - 1
      t2 <- ts_mean_test(a, lag)
      data.frame(forecast = nm[k], fz0 = mean(fz0_t[[k]][, i]), z1 = z1[1L], z1.p = 2 * stats::pnorm(-abs(z1[2L])), z2 = t2[["mean"]], z2.p = t2[["p.value"]])
    }))
  }
  dm <- NULL
  if (K > 1L) {
    rows <- list()
    for (k in 2:K) {
      add <- function(loss, d) rows[[length(rows) + 1L]] <<- data.frame(forecast = nm[k], loss = loss, t(ts_mean_test(d, lag)), stringsAsFactors = FALSE)
      add("aqs", qs_t[[1L]] - qs_t[[k]])
      if (!anyNA(u[[1L]]$logdens) && !anyNA(u[[k]]$logdens)) add("neglogscore", u[[k]]$logdens - u[[1L]]$logdens)
      for (i in seq_len(np)) if (!anyNA(fz0_t[[1L]][, i]) && !anyNA(fz0_t[[k]][, i])) add(paste0("fz0_", format(p[i])), fz0_t[[1L]][, i] - fz0_t[[k]][, i])
    }
    dm <- do.call(rbind, rows)
    names(dm)[names(dm) == "mean"] <- "diff"; names(dm)[names(dm) == "stat"] <- "statistic"
    rownames(dm) <- NULL
  }
  structure(list(summary = summ, profile = profile, var = var_tab, es = es_tab, dm = dm, levels = lev, weights = w, n = n, p = p),
            class = "tsscore")
}

#' @export
print.tsscore <- function(x, digits = 4, ...) {
  cat("Forecast scores for ", x$n, " observations\n\n", sep = "")
  s <- x$summary
  num <- vapply(s, is.numeric, TRUE)
  s[num] <- lapply(s[num], function(v) signif(v, digits))
  print(s, row.names = FALSE)
  cat("(all scores are losses: smaller is better; neglogscore = minus the mean log predictive density)\n")
  for (key in names(x$var)) {
    cat("\nVaR backtests, p = ", key, "\n", sep = "")
    v <- x$var[[key]]; nn <- vapply(v, is.numeric, TRUE); v[nn] <- lapply(v[nn], function(z) signif(z, digits))
    print(v, row.names = FALSE)
    e <- x$es[[key]]
    if (!all(is.na(e$fz0))) {
      cat("ES backtests, p = ", key, "\n", sep = "")
      nn <- vapply(e, is.numeric, TRUE); e[nn] <- lapply(e[nn], function(z) signif(z, digits))
      print(e, row.names = FALSE)
    }
  }
  if (!is.null(x$dm)) {
    cat("\nDiebold-Mariano tests against '", x$summary$forecast[1L], "' (positive: the other forecast is better)\n", sep = "")
    d <- x$dm; nn <- vapply(d, is.numeric, TRUE); d[nn] <- lapply(d[nn], function(z) signif(z, digits))
    print(d, row.names = FALSE)
  }
  invisible(x)
}

#' @param relative for the `plot` method: show the difference from the first forecast.
#' @rdname tsscore
#' @method plot tsscore
#' @export
plot.tsscore <- function(x, relative = FALSE, ...) {
  P <- x$profile
  if (relative) P <- P - P[, 1L]
  args <- list(...)
  K <- ncol(P)
  cols <- if (is.null(args$col)) seq_len(K) else rep_len(args$col, K)
  args$col <- NULL
  do.call(graphics::matplot, c(list(x$levels, P, type = "l", lty = 1, lwd = 2, col = cols,
    xlab = "level", ylab = if (relative) "quantile score minus first" else "mean quantile score"), args))
  if (relative) graphics::abline(h = 0, lty = 3)
  if (K > 1L) graphics::legend("top", colnames(P), col = cols, lwd = 2, bty = "n")
  invisible(x)
}
