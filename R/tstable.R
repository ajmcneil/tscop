# Tables of forecast scores: a tidy data frame from a tsscore object and a table of one loss for one or several score objects.

#' Scores of forecasts as a data frame
#'
#' Turns the result of [tsscore()] into a tidy data frame with one row for each forecast and loss.
#'
#' @details The losses are those in the `summary` of the score object: the average quantile score `aqs`, the negative
#' log score `neglogscore` and, if the forecasts contain risk measures, the Fissler-Ziegel losses `fz0_<p>`; a smaller loss
#' is better. The first forecast is the reference: the columns `change` and `change_pct` give the difference between the
#' score of a forecast and that of the reference and the percentage difference (`NA` if the score of the reference is not
#' positive, as for the negative log score, which can be negative); `statistic` and `p.value` belong to the
#' Diebold-Mariano test of the null hypothesis that the forecast is as accurate as the reference. As in [tsscore()],
#' a positive statistic (a negative `change`) means that the forecast is better than the reference. The rows of the
#' reference have `NA` in the test columns, as do tests that could not be computed (for example when a log predictive
#' density is infinite).
#'
#' @param x a `tsscore` object, the result of [tsscore()].
#' @param row.names,optional not used.
#' @param ... not used.
#' @return A data frame with the columns `forecast`, `loss`, `score`, `change`, `change_pct`, `statistic` and `p.value`,
#'   and the attribute `reference`, the name of the reference forecast.
#' @seealso [tsscore()], [tstable()]
#' @method as.data.frame tsscore
#' @export
#' @examples
#' set.seed(1)
#' m <- tsmod(marg("st", c(mu = 0, sigma = 1, df = 6)), armacop(ar = 0.5))
#' x <- tssim(m, 400)
#' g <- tsmod(marg("gauss", c(mu = 0, sigma = 1.3)), armacop(ar = 0.5))
#' good <- tsroll(m, x, start = 301, p = NULL)
#' bad <- tsroll(g, x, start = 301, p = NULL)
#' as.data.frame(tsscore(bad, good))
as.data.frame.tsscore <- function(x, row.names = NULL, optional = FALSE, ...) {
  s <- x$summary
  losses <- setdiff(names(s), "forecast")
  dm <- x$dm
  out <- do.call(rbind, lapply(losses, function(l) {
    score <- s[[l]]
    stat <- pv <- rep(NA_real_, length(score))
    d <- if (is.null(dm)) NULL else dm[dm$loss == l, , drop = FALSE]
    if (!is.null(d) && nrow(d) > 0L) {
      j <- match(d$forecast, s$forecast)
      stat[j] <- d$statistic
      pv[j] <- d$p.value
    }
    change <- score - score[1L]
    pct <- if (is.finite(score[1L]) && score[1L] > 0) 100 * change / score[1L] else rep(NA_real_, length(score))
    data.frame(forecast = s$forecast, loss = l, score = score, change = change, change_pct = pct,
               statistic = stat, p.value = pv, stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  attr(out, "reference") <- s$forecast[1L]
  out
}

#' Table of the scores of several forecasts
#'
#' Arranges one loss of the results of [tsscore()] in a table with a row for each forecast: its score, the change relative
#' to the reference forecast (the first one) and the p-value of the Diebold-Mariano test against the reference. The
#' results of several comparisons that use the same forecasts, for example with parameters held fixed or re-estimated,
#' can be put side by side.
#'
#' @details The argument `relative` decides how the change is expressed: as a percentage of the score of the reference
#' (`"percent"`, suitable for the positive losses `aqs` and `fz0_<p>`) or as a difference (`"difference"`, for the negative
#' log score, for which a percentage has no meaning). The default `"auto"` uses the percentage for the average quantile
#' score and the Fissler-Ziegel losses and the difference for the negative log score. The test is that of [tsscore()]:
#' it tests the difference between the loss of the forecast and that of the reference, and a small p-value means that the
#' forecasts differ in accuracy; a negative change means that the forecast is better. With several set-ups the forecasts
#' and the reference must be the same in all of them.
#'
#' @param x a `tsscore` object, or a named list of them; the names label the set-ups.
#' @param loss the loss to tabulate: `"aqs"` (default), `"neglogscore"` or, if the forecasts contain risk measures, a
#'   Fissler-Ziegel loss such as `"fz0_0.05"`.
#' @param relative how to express the change relative to the reference forecast; see Details.
#' @param scale a number by which the scores (and the differences) are multiplied, for example 1000 to avoid many
#'   decimals.
#' @return An object of class `tstable`, a data frame with the column `forecast` and, for each set-up, the columns
#'   `score`, `change` and `p.value` (with the name of the set-up as a suffix if there are several), and the attributes
#'   `loss`, `reference`, `relative` and `scale`. Printing rounds the numbers and describes the table in a header.
#' @seealso [tsscore()], [as.data.frame.tsscore()]
#' @export
#' @examples
#' set.seed(1)
#' m <- tsmod(marg("st", c(mu = 0, sigma = 1, df = 6)), armacop(ar = 0.5))
#' x <- tssim(m, 400)
#' g <- tsmod(marg("gauss", c(mu = 0, sigma = 1.3)), armacop(ar = 0.5))
#' good <- tsroll(m, x, start = 301, p = NULL)
#' bad <- tsroll(g, x, start = 301, p = NULL)
#' sc <- tsscore(bad, good)
#' tstable(sc)
#' tstable(sc, loss = "neglogscore")
#' tstable(list(first = sc, second = sc), scale = 1000)
tstable <- function(x, loss = "aqs", relative = c("auto", "percent", "difference"), scale = 1) {
  relative <- match.arg(relative)
  single <- inherits(x, "tsscore")
  setups <- if (single) list(x) else x
  if (!is.list(setups) || length(setups) == 0L || !all(vapply(setups, inherits, TRUE, "tsscore"))) {
    stop("'x' must be a tsscore object or a named list of tsscore objects", call. = FALSE)
  }
  if (!single && (is.null(names(setups)) || any(names(setups) == ""))) stop("the list of score objects must be named", call. = FALSE)
  if (!is.character(loss) || length(loss) != 1L) stop("'loss' must be the name of one loss", call. = FALSE)
  if (!is.numeric(scale) || length(scale) != 1L || !is.finite(scale)) stop("'scale' must be a number", call. = FALSE)
  tabs <- lapply(setups, function(sc) {
    d <- as.data.frame(sc)
    if (!(loss %in% d$loss)) stop("no loss '", loss, "' in the scores; available: ", paste(unique(d$loss), collapse = ", "), call. = FALSE)
    d[d$loss == loss, , drop = FALSE]
  })
  fc <- tabs[[1L]]$forecast
  if (!all(vapply(tabs, function(d) identical(d$forecast, fc), TRUE))) {
    stop("the forecasts, and in particular the reference, must be the same in all the score objects", call. = FALSE)
  }
  if (relative == "auto") relative <- if (loss == "aqs" || startsWith(loss, "fz0")) "percent" else "difference"
  out <- data.frame(forecast = fc, stringsAsFactors = FALSE)
  for (k in seq_along(tabs)) {
    d <- tabs[[k]]
    cols <- data.frame(score = scale * d$score,
                       change = if (relative == "percent") d$change_pct else scale * d$change,
                       p.value = d$p.value)
    if (!single) names(cols) <- paste0(names(cols), ".", names(tabs)[k])
    out <- cbind(out, cols)
  }
  structure(out, class = c("tstable", "data.frame"), loss = loss, reference = fc[1L], relative = relative, scale = scale)
}

#' @param digits number of significant digits when printing.
#' @param ... not used.
#' @rdname tstable
#' @method print tstable
#' @export
print.tstable <- function(x, digits = 4, ...) {
  scale <- attr(x, "scale")
  cat("Loss: ", attr(x, "loss"), if (scale != 1) paste0(" (scores multiplied by ", format(scale), ")"), "; smaller is better\n", sep = "")
  cat("change: ", if (attr(x, "relative") == "percent") "percentage change" else "difference", " relative to '", attr(x, "reference"), "'\n", sep = "")
  cat("p.value: Diebold-Mariano test against '", attr(x, "reference"), "'\n", sep = "")
  d <- as.data.frame(unclass(x), stringsAsFactors = FALSE)
  attributes(d)[c("loss", "reference", "relative", "scale")] <- NULL
  num <- vapply(d, is.numeric, TRUE)
  d[num] <- lapply(d[num], signif, digits = digits)
  print(d, row.names = FALSE)
  invisible(x)
}
