# Registry of the parametric margin families (see docs/marg-design-notes.md).
#
# A margin is a location-scale transform x = mu + sigma z of a standardised
# family whose shape parameters are collected in the named vector `s`. Every
# family supplies the log-density, distribution function (lower or upper tail)
# and quantile function (from a lower or upper tail probability) of the
# standardised variable z; `mg_*` helpers apply the location and scale.

# ---- symmetric standardised building blocks ---------------------------------

mg_gauss <- list(
  ld = function(z, s) stats::dnorm(z, log = TRUE),
  p = function(z, s, lt) stats::pnorm(z, lower.tail = lt),
  q = function(u, s, lt) stats::qnorm(u, lower.tail = lt)
)

# standard Laplace density exp(-|z|) / 2
mg_laplace <- list(
  ld = function(z, s) -abs(z) - log(2),
  p = function(z, s, lt) {
    if (!lt) z <- -z
    ifelse(z < 0, 0.5 * exp(z), 1 - 0.5 * exp(-z))
  },
  q = function(u, s, lt) {
    out <- ifelse(u < 0.5, log(2 * u), -log(2 * (1 - u)))
    if (lt) out else -out
  }
)

mg_st <- list(
  ld = function(z, s) stats::dt(z, s[["df"]], log = TRUE),
  p = function(z, s, lt) stats::pt(z, s[["df"]], lower.tail = lt),
  q = function(u, s, lt) stats::qt(u, s[["df"]], lower.tail = lt)
)

# Johnson SU: y = (x - mu) / sigma satisfies gamma + delta asinh(y) ~ N(0, 1),
# with delta > 0 (tail weight: smaller is heavier) and gamma real (skewness).
mg_jsu <- list(
  ld = function(y, s) {
    g <- s[["gamma"]]
    d <- s[["delta"]]
    log(d) - 0.5 * log1p(y^2) - 0.5 * log(2 * pi) - 0.5 * (g + d * asinh(y))^2
  },
  p = function(y, s, lt) stats::pnorm(s[["gamma"]] + s[["delta"]] * asinh(y), lower.tail = lt),
  q = function(u, s, lt) sinh((stats::qnorm(u, lower.tail = lt) - s[["gamma"]]) / s[["delta"]])
)

# Normal inverse Gaussian: z = skew * W + sqrt(W) * N with W inverse Gaussian of
# mean 1 and shape `shape`, so that, with a = shape and r = skew, the density of z is
#   sqrt(a (a + r^2)) / pi * exp(a + r z) K_1(sqrt((a + r^2)(a + z^2))) / sqrt(a + z^2)
# (the NIG density with alpha = sqrt(a + r^2), beta = r, delta = sqrt(a)).
# mu and sigma are then a pure location and scale; the mean of z is r and its
# variance 1 + r^2 / a. This is the ghyp alpha.bar parameterisation with the
# skewness gamma / sigma.
mg_nig_ld <- function(z, s) {
  a <- s[["shape"]]
  r <- s[["skew"]]
  x <- sqrt((a + r^2) * (a + z^2))
  0.5 * log(a + r^2) + 0.5 * log(a) - log(pi) + a + r * z - x +
    log(besselK(x, 1, expon.scaled = TRUE)) - 0.5 * log(a + z^2)
}
mg_nig <- list(
  ld = mg_nig_ld,
  p = function(z, s, lt) {
    mg_dens_cdf(function(t) exp(mg_nig_ld(t, s)), z, lt, s[["skew"]], min(1, sqrt(s[["shape"]])))
  },
  q = function(u, s, lt) {
    w0 <- min(1, sqrt(s[["shape"]]))
    dens <- function(t) exp(mg_nig_ld(t, s))
    mg_dens_quantile(function(t, l) mg_dens_cdf(dens, t, l, s[["skew"]], w0), dens, u, lt, s[["skew"]], w0)
  },
  r = function(n, s) {
    a <- s[["shape"]]
    # inverse Gaussian of mean 1 and shape a (Michael, Schucany and Haas, 1976)
    y <- stats::rnorm(n)^2
    w <- 1 + y / (2 * a) - sqrt(4 * a * y + y^2) / (2 * a)
    w <- ifelse(stats::runif(n) <= 1 / (1 + w), w, 1 / w)
    s[["skew"]] * w + sqrt(w) * stats::rnorm(n)
  }
)

# Fernandez-Steel skewing of a symmetric standardised family with skewing
# parameter gamma > 0 (gamma = 1 is the symmetric case; the mode is at 0):
#   f(z) = 2 / (gamma + 1 / gamma) * f0(gamma z) for z < 0, f0(z / gamma) for z >= 0.
# Both tails are computed from the corresponding tail of the base family, so
# that extreme upper-tail probabilities keep their accuracy.
mg_skew <- function(b) {
  list(
    ld = function(z, s) {
      g <- s[["gamma"]]
      log(2 / (g + 1 / g)) + ifelse(z < 0, b$ld(g * z, s), b$ld(z / g, s))
    },
    p = function(z, s, lt) {
      g <- s[["gamma"]]
      lo <- 2 / (g^2 + 1)
      hi <- 2 * g^2 / (g^2 + 1)
      neg <- z < 0
      out <- numeric(length(z))
      if (lt) {
        out[neg] <- lo * b$p(g * z[neg], s, TRUE)
        out[!neg] <- 1 - hi * b$p(z[!neg] / g, s, FALSE)
      } else {
        out[neg] <- 1 - lo * b$p(g * z[neg], s, TRUE)
        out[!neg] <- hi * b$p(z[!neg] / g, s, FALSE)
      }
      out
    },
    q = function(u, s, lt) {
      g <- s[["gamma"]]
      p0 <- 1 / (g^2 + 1) # P(Z < 0)
      out <- numeric(length(u))
      if (lt) {
        neg <- u < p0
        out[neg] <- b$q((g^2 + 1) * u[neg] / 2, s, TRUE) / g
        out[!neg] <- g * b$q((1 - u[!neg]) * (g^2 + 1) / (2 * g^2), s, FALSE)
      } else {
        pos <- u < 1 - p0
        out[pos] <- g * b$q(u[pos] * (g^2 + 1) / (2 * g^2), s, FALSE)
        out[!pos] <- b$q((1 - u[!pos]) * (g^2 + 1) / 2, s, TRUE) / g
      }
      out
    }
  )
}

# ---- the registry ----------------------------------------------------------
#
# types of the free-scale maps: "real" (identity), "pos" (log) and "df"
# (1 + exp, degrees of freedom above one). `start` gives starting values of the
# location and scale (and shapes) from the data; `parent` names the family
# whose fit starts the ladder, `from_parent` maps its parameters to the child's
# starting values.

mg_fam <- function(name, pars, types, base, start, parent = NA_character_, from_parent = NULL,
                   description = "") {
  names(types) <- pars
  list(name = name, pars = pars, types = types, base = base, start = start,
       parent = parent, from_parent = from_parent, description = description)
}

mg_robust <- function(x) {
  m <- stats::median(x)
  c(mu = m, mad = stats::mad(x, constant = 1.4826), sd = stats::sd(x), meanabs = mean(abs(x - m)))
}

.marg_families <- list(
  gauss = mg_fam("gauss", c("mu", "sigma"), c("real", "pos"), mg_gauss,
    start = function(x) c(mu = mean(x), sigma = stats::sd(x)),
    description = "Gaussian"),
  laplace = mg_fam("laplace", c("mu", "sigma"), c("real", "pos"), mg_laplace,
    start = function(x) {
      m <- stats::median(x)
      c(mu = m, sigma = mean(abs(x - m)))
    },
    description = "Laplace"),
  slaplace = mg_fam("slaplace", c("mu", "sigma", "gamma"), c("real", "pos", "pos"), mg_skew(mg_laplace),
    start = function(x) {
      m <- stats::median(x)
      c(mu = m, sigma = mean(abs(x - m)), gamma = 1)
    },
    parent = "laplace", from_parent = function(th) c(th, gamma = 1),
    description = "skew Laplace (Fernandez-Steel)"),
  st = mg_fam("st", c("mu", "sigma", "df"), c("real", "pos", "df"), mg_st,
    start = function(x) {
      m <- stats::median(x)
      df <- 8
      c(mu = m, sigma = stats::mad(x, constant = 1.4826) * 0.9, df = df)
    },
    parent = "gauss", from_parent = function(th) c(th, df = 8),
    description = "Student t"),
  jsu = mg_fam("jsu", c("mu", "sigma", "delta", "gamma"), c("real", "pos", "pos", "real"), mg_jsu,
    start = function(x) {
      d <- 2
      v <- 0.5 * (exp(1 / d^2) - 1) * (exp(1 / d^2) + 1) # variance of sinh(Z / d)
      c(mu = stats::median(x), sigma = stats::sd(x) / sqrt(v), delta = d, gamma = 0)
    },
    parent = "gauss",
    from_parent = function(th) {
      d <- 2
      v <- 0.5 * (exp(1 / d^2) - 1) * (exp(1 / d^2) + 1)
      c(mu = th[["mu"]], sigma = th[["sigma"]] / sqrt(v), delta = d, gamma = 0)
    },
    description = "Johnson SU"),
  nig = mg_fam("nig", c("mu", "sigma", "shape", "skew"), c("real", "pos", "pos", "real"), mg_nig,
    start = function(x) {
      m <- mean(x)
      sd <- stats::sd(x)
      zc <- (x - m) / sd
      ek <- max(mean(zc^4) - 3, 0.1)       # excess kurtosis 3 / shape for the symmetric NIG
      shape <- min(max(3 / ek, 0.1), 20)
      sk <- mean(zc^3)                     # skewness is approximately 3 skew / shape / sd(z)
      skew <- max(min(sk * shape / 3, 2), -2)
      c(mu = m - sd * skew / sqrt(1 + skew^2 / shape), sigma = sd / sqrt(1 + skew^2 / shape), shape = shape, skew = skew)
    },
    description = "normal inverse Gaussian"),
  sst = mg_fam("sst", c("mu", "sigma", "df", "gamma"), c("real", "pos", "df", "pos"),
    mg_skew(mg_st),
    start = function(x) {
      m <- stats::median(x)
      c(mu = m, sigma = stats::mad(x, constant = 1.4826) * 0.9, df = 8, gamma = 1)
    },
    parent = "st", from_parent = function(th) c(th, gamma = 1),
    description = "skew Student t (Fernandez-Steel)")
)

mg_family <- function(name) {
  fam <- .marg_families[[name]]
  if (is.null(fam)) {
    stop("unknown margin family '", name, "'; available: ", paste(names(.marg_families), collapse = ", "),
         call. = FALSE)
  }
  fam
}

# ---- free-scale maps ---------------------------------------------------------

mg_to_free <- function(th, types) {
  out <- numeric(length(th))
  for (i in seq_along(th)) {
    out[i] <- switch(types[[i]], real = th[[i]], pos = log(th[[i]]), df = log(max(th[[i]] - 1, 1e-3)))
  }
  out
}

mg_from_free <- function(z, types) {
  out <- numeric(length(z))
  for (i in seq_along(z)) {
    out[i] <- switch(types[[i]], real = z[[i]], pos = exp(z[[i]]), df = 1 + exp(z[[i]]))
  }
  out
}

# ---- functions of the original variable -----------------------------------------

mg_rand <- function(fam, th, n) {
  if (!is.null(fam$base$r)) {
    th[["mu"]] + th[["sigma"]] * fam$base$r(n, mg_shapes(fam, th))
  } else {
    mg_qf(fam, th, stats::runif(n))
  }
}

mg_shapes <- function(fam, th) th[setdiff(fam$pars, c("mu", "sigma"))]

mg_ldens <- function(fam, th, x) {
  z <- (x - th[["mu"]]) / th[["sigma"]]
  fam$base$ld(z, mg_shapes(fam, th)) - log(th[["sigma"]])
}

mg_cdf <- function(fam, th, x, lower.tail = TRUE) {
  z <- (x - th[["mu"]]) / th[["sigma"]]
  fam$base$p(z, mg_shapes(fam, th), lower.tail)
}

mg_qf <- function(fam, th, p, lower.tail = TRUE) {
  th[["mu"]] + th[["sigma"]] * fam$base$q(p, mg_shapes(fam, th), lower.tail)
}
