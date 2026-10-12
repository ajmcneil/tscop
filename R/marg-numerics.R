# Numerical distribution and quantile functions for margins that are defined
# through a density only (the generalised hyperbolic family and its special cases).
#
# The distribution function is obtained by Gauss-Legendre integration of the
# density over the intervals between the sorted evaluation points, cumulated
# from an anchor in the corresponding tail (a single call of `integrate` over an
# infinite range), so that both tails keep their relative accuracy: points at or
# below `center` accumulate the lower tail from the left, points above it the
# upper tail from the right. The quantile function interpolates the distribution
# function tabulated on a grid and polishes by Newton steps.

# Integrals of `dens` over the intervals [a_i, b_i] (a_i < b_i): composite
# 12-point Gauss-Legendre with panel widths proportional to the distance from
# the centre (but at least of the order w0), so that the narrow centre and the
# exponential tails are both resolved.
mg_panels <- function(dens, a, b, center, w0, nodes = 12L) {
  if (length(a) == 0L) return(numeric(0))
  g <- gl_nodes(nodes)
  dist <- pmin(abs(a - center), abs(b - center))
  cross <- (a - center) * (b - center) <= 0
  dist[cross] <- 0
  h <- 0.5 * pmax(w0, dist)
  nsub <- pmax(1, ceiling((b - a) / h))
  idx <- rep(seq_along(a), nsub)
  k <- sequence(nsub) - 1
  wd <- ((b - a) / nsub)[idx]
  left <- a[idx] + k * wd
  x <- outer(wd, g$s) + left
  v <- matrix(dens(as.vector(x)), nrow = length(wd))
  v[!is.finite(v)] <- 0
  as.vector(rowsum(wd * as.vector(v %*% g$w), idx))
}

# Distribution function of the density `dens` at z: lower tail (lt = TRUE) or
# upper tail. `center` splits the points between the two accumulations; `w0` is
# the width of the central region.
mg_dens_cdf <- function(dens, z, lt, center, w0) {
  out <- rep(NA_real_, length(z))
  fin <- is.finite(z)
  out[!is.na(z) & z == -Inf] <- if (lt) 0 else 1
  out[!is.na(z) & z == Inf] <- if (lt) 1 else 0
  if (!any(fin)) return(out)
  ux <- sort(unique(z[fin]))
  nu <- length(ux)
  low <- ux <= center
  F <- numeric(nu) # lower-tail probability for the low points, upper-tail for the high ones
  anchor <- function(f, lower, upper) {
    r <- tryCatch(stats::integrate(f, lower, upper, rel.tol = 1e-11, abs.tol = 0, subdivisions = 500L),
                  error = function(e) NULL)
    if (is.null(r)) 0 else max(r$value, 0)
  }
  if (any(low)) {
    il <- which(low)
    F[il[1L]] <- anchor(dens, -Inf, ux[il[1L]])
    if (length(il) > 1L) {
      pn <- mg_panels(dens, ux[il[-length(il)]], ux[il[-1L]], center, w0)
      F[il[-1L]] <- F[il[1L]] + cumsum(pn)
    }
  }
  if (any(!low)) {
    ih <- which(!low)
    top <- ih[length(ih)]
    F[top] <- anchor(dens, ux[top], Inf)
    if (length(ih) > 1L) {
      pn <- mg_panels(dens, ux[ih[-length(ih)]], ux[ih[-1L]], center, w0)
      F[rev(ih)[-1L]] <- F[top] + cumsum(rev(pn))
    }
  }
  res <- if (lt) ifelse(low, F, 1 - F) else ifelse(low, 1 - F, F)
  out[fin] <- res[match(z[fin], ux)]
  pmin(pmax(out, 0), 1)
}

# Quantile function: solve P(Z <= z) = u (lt = TRUE) or P(Z > z) = u (lt = FALSE).
# `pdf(z, lt)` and `dens(z)` are the distribution function (either tail) and
# density.
mg_dens_quantile <- function(pdf, dens, u, lt, center, w0) {
  out <- rep(NA_real_, length(u))
  ok <- !is.na(u) & u > 0 & u < 1
  out[!is.na(u) & u == 0] <- if (lt) -Inf else Inf
  out[!is.na(u) & u == 1] <- if (lt) Inf else -Inf
  if (!any(ok)) return(out)
  # lower-tail target probability and a flag for solving in the upper tail
  up <- if (lt) u > 0.5 else u < 0.5       # solve on the upper tail, with target su
  su <- if (lt) 1 - u else u
  pl <- if (lt) u else 1 - u
  tgt <- ifelse(up, su, pl)
  utgt <- u[ok]
  tg <- tgt[ok]
  upk <- up[ok]
  # bracket the extreme targets by doubling the distance from the centre
  tmin <- min(tg)
  dl <- max(w0, 1)
  zL <- center - dl
  for (i in 1:60) {
    if (pdf(zL, TRUE) <= 0.5 * tmin) break
    dl <- 2 * dl
    zL <- center - dl
  }
  dr <- max(w0, 1)
  zR <- center + dr
  for (i in 1:60) {
    if (pdf(zR, FALSE) <= 0.5 * tmin) break
    dr <- 2 * dr
    zR <- center + dr
  }
  grid <- sort(unique(c(center + sinh(seq(asinh(zL - center), asinh(zR - center), length.out = 500L)))))
  Fg <- mg_dens_cdf(dens, grid, TRUE, center, w0)
  Sg <- mg_dens_cdf(dens, grid, FALSE, center, w0)
  z <- numeric(length(tg))
  for (side in c(FALSE, TRUE)) {
    sel <- upk == side
    if (!any(sel)) next
    if (!side) { # lower tail: F increasing in the grid
      lf <- log(pmax(Fg, 1e-300))
      z[sel] <- stats::approx(lf, grid, xout = log(tg[sel]), rule = 2, ties = "ordered")$y
    } else {     # upper tail: S decreasing in the grid
      ls <- log(pmax(Sg, 1e-300))
      z[sel] <- stats::approx(rev(ls), rev(grid), xout = log(tg[sel]), rule = 2, ties = "ordered")$y
    }
  }
  # Newton polish on the tail in which each point was solved
  for (it in 1:6) {
    Fz <- numeric(length(z))
    if (any(!upk)) Fz[!upk] <- mg_dens_cdf(dens, z[!upk], TRUE, center, w0)
    if (any(upk)) Fz[upk] <- mg_dens_cdf(dens, z[upk], FALSE, center, w0)
    fz <- dens(z)
    step <- (Fz - tg) / fz
    step[!is.finite(step)] <- 0
    step[upk] <- -step[upk]
    # damp: never move by more than half the distance to the centre plus a unit
    lim <- 0.5 * abs(z - center) + 1
    step <- pmax(pmin(step, lim), -lim)
    z <- z - step
    if (max(abs(step) / (1 + abs(z))) < 1e-13) break
  }
  out[ok] <- z
  out
}
