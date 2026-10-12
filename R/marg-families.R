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
