# Pair copula specifications for tied d-vine models.

.pc_families <- c("gauss", "t", "gumbel", "clayton", "joe", "frank", "bb1", "ast", "ivt")
# families without rotation arguments: radially symmetric ones, and the IVT
# pair copula (positive dependence only, the base copula is fixed)
.pc_radial <- c("gauss", "t", "frank", "ivt")
# base copulas of an inverse-v-transformed pair copula: positive dependence and
# a finite density at the origin ("clayton180" is the survival Clayton)
.pc_ivt_bases <- c("ast", "joe", "gumbel", "clayton180", "frank", "gauss")
.pc_aux_families <- c("t", "bb1")
# families that tend to independence as Kendall's tau -> 0 with auxiliary
# parameters fixed, hence admissible as the tail family
.pc_tail_families <- c("gauss", "frank", "clayton", "gumbel", "joe", "ast", "ivt")

#' Pair copula specifications
#'
#' Specification of the pair copula used at a lag of a tied [sdvinecop-class]
#' model. The copula is calibrated at each lag to the Kendall's tau value of
#' the model's Kendall partial autocorrelation function, so the specification
#' needs only the family, the rotations to use under positive and negative
#' dependence, and (for two-parameter families) the auxiliary parameter.
#'
#' @slot family one of `"gauss"`, `"t"`, `"gumbel"`, `"clayton"`, `"joe"`,
#'   `"frank"`, `"bb1"`, `"ast"` and `"ivt"`. `"ivt"` is the inverse-v-transformed
#'   pair copula of Dias, Han and McNeil (2027), `c(u, v) = c*(V1(u), V2(v))`
#'   with linear v-transforms whose fulcrums are parameters of the d-vine model
#'   (`vt` in [sdvinecop-class]) and a base copula with positive dependence
#'   given by `base`; `"ast"` is (the absolute spherical t copula of Dias,
#'   Han and McNeil, from the udp package: one parameter, positive dependence
#'   and upper tail dependence, like Gumbel and Joe; calibrated to Kendall's
#'   tau through its degrees-of-freedom parameter; roughly ten times slower to
#'   evaluate than the `rvinecopulib` families).
#' @slot pos rotation (0 or 180) used when Kendall's tau is non-negative.
#' @slot neg rotation (90 or 270) used when Kendall's tau is negative.
#' @slot base for `"ivt"`: the base copula, one of `"ast"`, `"joe"`, `"gumbel"`,
#'   `"clayton180"` (survival Clayton), `"frank"` (no tail dependence, cheaper
#'   than the tail-dependent bases) and `"gauss"` (Gaussian, the cheapest, a
#'   base for weakly dependent high lags); `NA` otherwise.
#' @slot aux auxiliary parameter: degrees of freedom for `"t"` (at least 2)
#'   or the second parameter (`delta` in rvinecopulib) for `"bb1"` (at least 1); `NA` for other families, and `NA`
#'   for `"t"`/`"bb1"` if it has not been given yet (a starting value is then
#'   chosen when the model is fitted).
#'
#' @references Dias, A., Han, J. and McNeil, A. J. (2027). GARCH copulas, v-transforms and
#'   D-vines for stochastic volatility. *Journal of Multivariate Analysis*, 217,
#'   105695. \doi{10.1016/j.jmva.2026.105695}
#'
#' @exportClass pc
setClass("pc", slots = c(family = "character", pos = "numeric",
                         neg = "numeric", aux = "numeric", base = "character"))

setValidity("pc", function(object) {
  if (length(object@family) != 1L || !(object@family %in% .pc_families)) {
    return(paste0("'family' must be one of ", paste(.pc_families, collapse = ", ")))
  }
  if (length(object@pos) != 1L || length(object@neg) != 1L || length(object@aux) != 1L) {
    return("'pos', 'neg' and 'aux' must be single numbers")
  }
  if (object@family %in% .pc_radial) {
    if (object@pos != 0 || object@neg != 0) return("radially symmetric families have pos = neg = 0")
  } else {
    if (!(object@pos %in% c(0, 180))) return("'pos' must be 0 or 180")
    if (!(object@neg %in% c(90, 270))) return("'neg' must be 90 or 270")
  }
  if (length(object@base) != 1L) return("'base' must be a single string")
  if (object@family == "ivt") {
    if (is.na(object@base) || !(object@base %in% .pc_ivt_bases)) {
      return(paste0("an \"ivt\" pair copula needs a base copula, one of ", paste(.pc_ivt_bases, collapse = ", ")))
    }
  } else if (!is.na(object@base)) {
    return("only the \"ivt\" family has a base copula")
  }
  if (object@family %in% .pc_aux_families) {
    if (!is.na(object@aux)) {
      if (object@family == "t" && object@aux < 2) return("degrees of freedom must be at least 2")
      if (object@family == "bb1" && object@aux < 1) return("the second BB1 parameter 'aux' must be at least 1")
    }
  } else if (!is.na(object@aux)) {
    return(paste0("family '", object@family, "' has no auxiliary parameter"))
  }
  TRUE
})

#' Constructor for pair copula specifications
#'
#' @param family family name (case-insensitive): `"gauss"`, `"t"`,
#'   `"gumbel"`, `"clayton"`, `"joe"`, `"frank"`, `"bb1"` or `"ast"`.
#' @param pos rotation (0 or 180) used when Kendall's tau at the lag is
#'   non-negative. Ignored for radially symmetric families (Gauss, t, Frank).
#' @param neg rotation (90 or 270) used when Kendall's tau at the lag is
#'   negative. Ignored for radially symmetric families.
#' @param aux auxiliary parameter: degrees of freedom for `"t"`, the second BB1 parameter (`delta` in rvinecopulib) for
#'   `"bb1"`; `NULL` for the other families.
#' @param base for `family = "ivt"` the base copula: `"ast"`, `"joe"`,
#'   `"gumbel"`, `"clayton180"` (the survival Clayton copula), `"frank"` or `"gauss"`
#'   (the Gaussian copula, whose inverse-v-transformed version has a
#'   cross-shaped density; cheap to evaluate, so useful for the weakly
#'   dependent high lags); `NULL` for the other families.
#'
#' @details Both rotations are always stored. Which one is used at a lag is
#'   decided by the sign of Kendall's tau at the current parameter values, so
#'   the user does not need to know the sign in advance. The rotation
#'   convention is that of [rvinecopulib::bicop_dist()] applied to the pair
#'   (earlier variable, later variable).
#'
#' @return An object of class [pc-class].
#' @references Dias, A., Han, J. and McNeil, A. J. (2027). GARCH copulas, v-transforms and
#'   D-vines for stochastic volatility. *Journal of Multivariate Analysis*, 217,
#'   105695. \doi{10.1016/j.jmva.2026.105695}
#'
#' @export
#'
#' @examples
#' pc("gumbel")
#' pc("clayton", pos = 180, neg = 270)
#' pc("t", aux = 6)
#' pc("ivt", base = "ast")
#' rep(pc("ivt", base = "ast"), 3)   # the same specification at three lags
#' c(rep(pc("ivt", base = "ast"), 3), rep(pc("gauss"), 2))
pc <- function(family, pos = 0, neg = 90, aux = NULL, base = NULL) {
  if (!is.character(family) || length(family) != 1L) stop("'family' must be a single name")
  family <- tolower(family)
  family <- switch(family, gaussian = "gauss", normal = "gauss", student = "t",
                   student_t = "t", family)
  if (!(family %in% .pc_families)) {
    stop("'family' must be one of ", paste(.pc_families, collapse = ", "))
  }
  if (family %in% .pc_radial) {
    pos <- 0
    neg <- 0
  }
  if (is.null(aux)) aux <- NA_real_
  if (is.null(base)) {
    base <- NA_character_
  } else {
    base <- tolower(base)
  }
  if (!is.numeric(pos) || !is.numeric(neg) || !is.numeric(aux)) stop("'pos', 'neg' and 'aux' must be numeric")
  new("pc", family = family, pos = as.numeric(pos), neg = as.numeric(neg),
      aux = as.numeric(aux), base = as.character(base))
}

# Coerce a family name or pc object to a pc object.
as_pc <- function(x) {
  if (is(x, "pc")) return(x)
  if (is.character(x) && length(x) == 1L) return(pc(x))
  stop("expected a pc() specification or a family name")
}

setClassUnion("pcOrNULL", c("pc", "NULL"))
setClassUnion("armacopOrNULL", c("armacop", "NULL"))

#' @rdname tscop-methods
#' @export
setMethod("show", "pc", function(object) {
  rot <- if (object@family == "ivt") paste0(" with base ", object@base)
         else if (object@family %in% .pc_radial) ""
         else sprintf(" (positive: %g, negative: %g)", object@pos, object@neg)
  aux <- if (is.na(object@aux)) "" else sprintf(", aux = %g", object@aux)
  cat("pair copula specification: ", object@family, rot, aux, "\n", sep = "")
  invisible(object)
})

#' @rdname pc
#' @param x a [pc-class] specification, for `rep`.
#' @param ... arguments of [base::rep()] (`times`, `each`, `length.out`) for `rep`.
#' @details `rep(x, ...)` repeats a specification and returns a list of
#'   [pc-class] objects, the form taken by the `pairs` argument of [sdvinecop()].
#' @export
setMethod("rep", "pc", function(x, ...) base::rep(list(x), ...))
