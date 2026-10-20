# Plots of the pair copulas of a d-vine copula process: contours of the density of the pair copula at selected lags,
# optionally with the data pairs of that lag (the pairs the pair copula describes, conditioned on the observations
# in between through the h-functions) and, for inverse-v-transformed pair copulas, the two v-transforms in marginal strips.

#' Plot the pair copulas of a d-vine copula process
#'
#' Draws the density of the pair copula \eqn{c_k(u, v)} at the chosen lags \eqn{k} of a [sdvinecop-class] model as filled
#' contour plots in the style of the plot method for `bsicopula` objects of the udp package. For lag \eqn{k} this is the
#' copula of \eqn{(U_t, U_{t+k})} given the intermediate observations: the partial copula. For models with inverse-v-transformed
#' (`"ivt"`) pair copulas the contours are cross-shaped, with their centre at the fulcrums of the two v-transforms.
#'
#' Contour levels are by default quantiles of the density values pooled over all the panels, so that the panels are comparable;
#' the filled background is coloured by the rank of the density, as in udp, and is therefore an ordinal scale. Pass `levels` to
#' use fixed density values instead (for example to compare two models).
#'
#' With `data` the pairs of the data at lag \eqn{k} can be overlaid: at lag 1 these are the consecutive observations, at higher lags
#' the pairs conditioned on the observations in between through the h-functions of the pair copulas of the lower lags. They
#' are the pairs that the pair copula of that lag describes. With `strips = TRUE`, panels of inverse-v-transformed pair copulas get
#' two marginal strips showing the v-transforms \eqn{V_1(u)} (above) and \eqn{V_2(v)} (to the right) and the fulcrums are
#' marked in the panel.
#'
#' @param x a [sdvinecop-class] model with specified coefficients, a fitted model ([sdvinecopfit-class], whose copula data
#'   are used for `data` unless given) or a fitted [tsmod-class] model with a `sdvinecop` copula (the data are the probability
#'   integral transforms of the margin).
#' @param y `NULL`/missing, or a second [sdvinecop-class] model (or fitted model) to draw next to `x` for the same lags: the panels
#'   are then arranged with one row per lag and one column per model, with common contour levels, for example to compare a
#'   linear and a nonlinear v-transform.
#' @param lags the lags to draw, at most the effective lag `K` (default `1:4`).
#' @param data copula data in (0, 1) to overlay as points (default none; for fitted models `TRUE` uses the data of the fit).
#' @param points the lags at which the data points are drawn (default all of `lags` when `data` is given).
#' @param strips logical; draw the v-transforms of inverse-v-transformed pair copulas in marginal strips.
#' @param n number of grid points per axis.
#' @param fill logical; colour the background by the rank of the density (default). Use `FALSE` for contour lines only, which
#'   shows overlaid data points better.
#' @param probs probabilities of the density quantiles used as contour levels when `levels` is `NULL`.
#' @param levels contour levels (density values) instead of `probs`.
#' @param col colours of the filled background, as in the `bsicopula` method.
#' @param pt.col colour of the data points.
#' @param labels names of the two models when `y` is given.
#' @param main titles of the panels (default `"lag k"` with the Kendall's tau of the pair copula).
#' @param ... further arguments to [graphics::image()].
#' @return Invisibly, a list with the contour `levels`, the `lags` and the Kendall's tau `tau` of the pair copulas at the lags.
#' @seealso [ktauhat()], [sdvinecop()], [pc()]
#' @rdname plot-sdvinecop
#' @importFrom grDevices hcl.colors adjustcolor
#' @export
#' @examples
#' kp <- arfactcop(lambda = 0.95, w = 0.3)
#' m <- sdvinecop(kp, tail = pc("ivt", base = "joe"), vt = c(0.45, 0.65), vttype = "v2p",
#'                vtshape = c(0.7, 1.4), maxlag = 12, tautol = 1e-3)
#' plot(m)                                  # lags 1 to 4
#' plot(m, lags = 1, strips = TRUE)         # with the two v-transforms in strips
#' mg <- marg("st", c(mu = 0, sigma = 1, df = 6))
#' u <- mdf(mg, tssim(tsmod(mg, m), 1000))
#' plot(m, lags = 1, data = u, strips = TRUE)
#' mlin <- sdvinecop(kp, tail = pc("ivt", base = "joe"), vt = c(0.3, 0.75), maxlag = 12, tautol = 1e-3)
#' plot(mlin, m, lags = 1, data = u, strips = TRUE, labels = c("linear", "v2p"))   # side by side
setMethod("plot", c(x = "sdvinecop", y = "missing"), function(x, y, lags = 1:4, data = NULL, points = NULL, strips = FALSE, n = 150, fill = TRUE,
         probs = c(seq(0.5, 0.9, 0.1), 0.95, 0.99, 0.999), levels = NULL,
         col = grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
         pt.col = grDevices::adjustcolor("grey25", 0.35), main = NULL, labels = NULL, ...) {
  dv_plot_pairs(list(x), lags, data, points, strips, n, fill, probs, levels, col, pt.col, main, labels, ...)
})

#' @rdname plot-sdvinecop
#' @export
setMethod("plot", c(x = "sdvinecop", y = "sdvinecop"), function(x, y, lags = 1:4, data = NULL, points = NULL, strips = FALSE, n = 150, fill = TRUE,
         probs = c(seq(0.5, 0.9, 0.1), 0.95, 0.99, 0.999), levels = NULL,
         col = grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
         pt.col = grDevices::adjustcolor("grey25", 0.35), main = NULL, labels = NULL, ...) {
  dv_plot_pairs(list(x, y), lags, data, points, strips, n, fill, probs, levels, col, pt.col, main, labels, ...)
})

#' @rdname plot-sdvinecop
#' @export
setMethod("plot", c(x = "sdvinecopfit", y = "missing"), function(x, y, lags = 1:4, data = NULL, points = NULL, strips = FALSE, n = 150, fill = TRUE,
         probs = c(seq(0.5, 0.9, 0.1), 0.95, 0.99, 0.999), levels = NULL,
         col = grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
         pt.col = grDevices::adjustcolor("grey25", 0.35), main = NULL, labels = NULL, ...) {
  if (isTRUE(data)) data <- x@data
  dv_plot_pairs(list(x@model), lags, data, points, strips, n, fill, probs, levels, col, pt.col, main, labels, ...)
})

#' @rdname plot-sdvinecop
#' @export
setMethod("plot", c(x = "sdvinecopfit", y = "sdvinecopfit"), function(x, y, lags = 1:4, data = NULL, points = NULL, strips = FALSE, n = 150, fill = TRUE,
         probs = c(seq(0.5, 0.9, 0.1), 0.95, 0.99, 0.999), levels = NULL,
         col = grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
         pt.col = grDevices::adjustcolor("grey25", 0.35), main = NULL, labels = NULL, ...) {
  if (isTRUE(data)) data <- x@data
  dv_plot_pairs(list(x@model, y@model), lags, data, points, strips, n, fill, probs, levels, col, pt.col, main, labels, ...)
})

#' @rdname plot-sdvinecop
#' @export
setMethod("plot", c(x = "tsmodfit", y = "missing"), function(x, y, lags = 1:4, data = NULL, points = NULL, strips = FALSE, n = 150, fill = TRUE,
         probs = c(seq(0.5, 0.9, 0.1), 0.95, 0.99, 0.999), levels = NULL,
         col = grDevices::hcl.colors(100, "YlOrRd", rev = TRUE),
         pt.col = grDevices::adjustcolor("grey25", 0.35), main = NULL, labels = NULL, ...) {
  cm <- if (is(x@copfit, "sdvinecopfit")) x@copfit@model else x@copfit
  if (!is(cm, "sdvinecop")) stop("the copula of the model is not a d-vine (sdvinecop) copula", call. = FALSE)
  if (isTRUE(data)) data <- tm_pits(x@margfit, x@data)
  dv_plot_pairs(list(cm), lags, data, points, strips, n, fill, probs, levels, col, pt.col, main, labels, ...)
})

# The v-transform of the variable `w` (1 = earlier, 2 = later) of an inverse-v-transformed pair copula, NULL for other pairs.
dv_pair_vt <- function(pr, w) {
  if (dv_is_ivtnl(pr)) function(x) dv_v2p(x, pr$delta[w], pr$kappa[w])
  else if (dv_is_ivt(pr)) function(x) dv_vlin(x, pr$delta[w])
  else NULL
}

# `models`: a list of one or two sdvinecop models. One model: the panels (one per lag) are arranged on a grid; two models: one row per
# lag and one column per model.
dv_plot_pairs <- function(models, lags, data, points, strips, n, fill, probs, levels, col, pt.col, main, labels, ...) {
  nm <- length(models)
  pairs <- lapply(models, function(m) {
    pr <- dv_resolve(m)
    if (is.null(pr)) stop("the model is infeasible at its current parameter values", call. = FALSE)
    pr
  })
  K <- min(vapply(pairs, length, 1L))
  if (!is.numeric(lags) || length(lags) < 1L || anyNA(lags) || any(lags != round(lags)) || any(lags < 1L) || any(lags > K)) {
    stop(sprintf("'lags' must be whole numbers between 1 and the effective lag K = %d", K), call. = FALSE)
  }
  lags <- as.integer(lags)
  np <- length(lags)
  g <- (seq_len(n) - 0.5) / n
  G <- expand.grid(u = g, v = g)
  dens <- lapply(seq_len(nm), function(j) lapply(lags, function(k) matrix(dv_dens(pairs[[j]][[k]], G$u, G$v), n, n)))
  if (is.null(levels)) levels <- sort(unique(stats::quantile(unlist(dens), probs, names = FALSE)))
  if (is.null(points)) points <- if (is.null(data)) integer(0) else lags
  pts <- vector("list", nm)
  if (length(points) > 0L) {
    if (is.null(data)) stop("'data' is needed to draw points", call. = FALSE)
    if (!is.numeric(data) || anyNA(data) || any(data <= 0 | data >= 1)) stop("'data' must be copula data strictly inside (0, 1)", call. = FALSE)
    if (!all(points %in% lags)) stop("'points' must be a subset of 'lags'", call. = FALSE)
    kmax <- max(points)
    if (length(data) < kmax + 2L) stop("too few data for the requested lags", call. = FALSE)
    for (j in seq_len(nm)) pts[[j]] <- dv_level_tau_r(pairs[[j]][seq_len(kmax)], as.numeric(data), keep = kmax)$pairs
  }
  tau <- vapply(seq_len(nm), function(j) vapply(lags, function(k) tryCatch(dv_pair_ktau(pairs[[j]][[k]]), error = function(e) NA_real_), 0), numeric(np))
  tau <- matrix(tau, np, nm)
  if (is.null(labels)) labels <- if (nm > 1L) paste("model", seq_len(nm)) else ""
  # panels in row-major order of the grid: (lag i, model j)
  grid <- if (nm == 1L) { nc <- if (np <= 3L) np else ceiling(sqrt(np)); cbind(i = seq_len(np), j = 1L, r = (seq_len(np) - 1L) %/% nc, c = (seq_len(np) - 1L) %% nc) }
          else { nc <- nm; cbind(i = rep(seq_len(np), each = nm), j = rep(seq_len(nm), np), r = rep(seq_len(np) - 1L, each = nm), c = rep(seq_len(nm) - 1L, np)) }
  nr <- max(grid[, "r"]) + 1L; nc <- max(grid[, "c"]) + 1L
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op))
  if (strips) {
    ids <- matrix(0L, 2L * nr, 2L * nc)
    for (p in seq_len(nrow(grid))) {
      r <- grid[p, "r"]; cc <- grid[p, "c"]
      ids[2L * r + 1L, 2L * cc + 1L] <- 3L * p - 2L
      ids[2L * r + 2L, 2L * cc + 1L] <- 3L * p - 1L
      ids[2L * r + 2L, 2L * cc + 2L] <- 3L * p
    }
    graphics::layout(ids, widths = rep(c(4, 1.1), nc), heights = rep(c(1.1, 4), nr))
  } else {
    graphics::par(mfrow = c(nr, nc))
  }
  xs <- seq(0, 1, length.out = 400)
  cv <- c("dodgerblue3", "darkorange2")
  dots <- list(...)
  for (p in seq_len(nrow(grid))) {
    i <- grid[p, "i"]; j <- grid[p, "j"]; k <- lags[i]; pr <- pairs[[j]][[k]]
    v1 <- dv_pair_vt(pr, 1L); v2 <- dv_pair_vt(pr, 2L)
    use_strip <- strips && !is.null(v1)
    ttl <- if (!is.null(main)) rep_len(main, nrow(grid))[p]
           else sprintf("lag %d%s  (tau = %.3f)", k, if (nm > 1L) paste0(", ", labels[j]) else "", tau[i, j])
    if (strips) {   # the strip panels are always created so that the layout stays regular
      graphics::par(mar = c(0.4, 4, 2.5, 0.5))
      if (use_strip) {
        graphics::plot(xs, v1(xs), type = "l", col = cv[1L], lwd = 2, xaxt = "n", yaxt = "n", xlim = c(0, 1), ylim = c(0, 1), xaxs = "i",
                       ylab = "V1(u)", main = ttl)
        graphics::axis(2, at = c(0, 1))
      } else graphics::plot.new()
    }
    graphics::par(mar = if (strips) c(4, 4, 0.4, 0.5) else c(4, 4, 2.5, 1))
    z <- matrix(rank(pmax(dens[[j]][[i]], min(levels))) / length(dens[[j]][[i]]), n, n)
    if (fill) {
      do.call(graphics::image, c(list(x = g, y = g, z = z, col = col, xlab = "u (earlier)", ylab = "v (later)", xaxs = "i", yaxs = "i",
                                      main = if (strips) "" else ttl), dots))
    } else {
      do.call(graphics::plot, c(list(NA, xlim = c(0, 1), ylim = c(0, 1), xlab = "u (earlier)", ylab = "v (later)", xaxs = "i", yaxs = "i",
                                     main = if (strips) "" else ttl), dots))
    }
    graphics::contour(g, g, dens[[j]][[i]], levels = levels, drawlabels = FALSE, add = TRUE, col = if (fill) "grey20" else "firebrick")
    if (k %in% points) graphics::points(pts[[j]][[k]][, 1L], pts[[j]][[k]][, 2L], pch = 16, cex = 0.35, col = pt.col)
    if (use_strip) {
      graphics::abline(v = pr$delta[1L], lty = 3, col = cv[1L]); graphics::abline(h = pr$delta[2L], lty = 3, col = cv[2L])
      graphics::par(mar = c(4, 1.2, 0.4, 3.5))
      graphics::plot(v2(xs), xs, type = "l", col = cv[2L], lwd = 2, xaxt = "n", yaxt = "n", xlim = c(0, 1), ylim = c(0, 1), yaxs = "i", xlab = "V2(v)")
      graphics::axis(1, at = c(0, 1)); graphics::axis(4)
    } else if (strips) graphics::plot.new()
  }
  invisible(list(levels = levels, lags = lags, tau = if (nm == 1L) drop(tau) else tau))
}
