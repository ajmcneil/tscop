#' Fit a time series copula model
#'
#' Generic for maximum likelihood estimation of the time series copula models
#' in tscop. Methods exist for [armacop-class] objects.
#'
#' @param U a numeric vector or time series of copula data in (0, 1):
#'   pseudo-observations (ranks), or probability integral transforms from
#'   fitted marginal models (set `pseudo = FALSE`). `NA` values are allowed.
#'   Values equal to 0 or 1 are rejected; see [udp::boundaryadjust()].
#' @param object the model to fit: a model structure, or a known model whose
#'   coefficients are used as starting values.
#' @param ... further arguments passed to methods.
#'
#' @return A fitted model object.
#' @export
setGeneric("tsfit", function(U, object, ...) standardGeneric("tsfit"))

#' Simulate from a time series copula model
#'
#' Generic for simulation from the time series copula models in tscop.
#' Methods exist for [armacop-class] and [armacopfit-class] objects.
#'
#' @param object a model object.
#' @param n length of the simulated series.
#' @param ... further arguments passed to methods.
#'
#' @return A numeric vector of copula data in (0, 1) of length `n`.
#' @export
setGeneric("tssim", function(object, n, ...) standardGeneric("tssim"),
           signature = "object")

#' Kendall's tau of the pair copulas of a time series copula model
#'
#' Generic returning Kendall's tau of the pair copulas at lags
#' `1, ..., lagmax` in the D-vine representation of the process. For an
#' [armacop-class] or [arfactcop-class] model the pair copulas are Gaussian
#' with parameters equal to the partial autocorrelations of the process, so
#' \eqn{\tau_k = (2/\pi)\arcsin(\alpha_k)} with \eqn{\alpha_k} the lag-`k`
#' partial autocorrelation.
#'
#' @param object a model object.
#' @param lagmax maximum lag.
#' @param ... further arguments passed to methods.
#'
#' @return A numeric vector of length `lagmax`.
#' @export
setGeneric("ktau", function(object, lagmax = 20, ...) standardGeneric("ktau"),
           signature = "object")

#' Residuals of time series copula models
#'
#' Residuals of a known model or of a fitted model object, iid uniform (or
#' standard normal on the `"normal"` scale) if the model is correct. If the
#' data are a `ts`, `zoo` or `xts` object the time index is restored.
#' \itemize{
#'   \item [armacop-class] and [arfactcop-class]: the standardised one-step
#'   prediction errors of the normal scores \eqn{\Phi^{-1}(u_t)}, on the
#'   default `"uniform"` scale their `pnorm`. Missing observations give `NA`.
#'   \item [sdvinecop-class]: the forward Rosenblatt transforms
#'   \eqn{R_{\min(t-1,K)}(u_t; u_{t-1}, \ldots)}, with lower-order
#'   Rosenblatt functions for the first \eqn{K} values. The data must be
#'   complete.
#'   \item [vtscop-class]: the conditional probability integral transforms
#'   \eqn{F(u_t \mid u_{t-1}, \ldots)}. The conditional distribution of
#'   \eqn{U_t} is that of the base process at \eqn{V(U_t)} passed through the
#'   stochastic inverse of the v-transform: closed form for linear
#'   v-transforms, otherwise a one-dimensional Gauss-Legendre integral. For a
#'   d-vine base with a nonlinear v-transform the residuals are computed one
#'   time point at a time and can be slow.
#'   \item [marg-class]: the probability integral transforms
#'   \eqn{F(x_t)} of the observations (uniform if the margin is correct; the
#'   `qnorm` of them on the `"normal"` scale).
#'   \item [tsmod-class]: the residuals of the copula process applied to the
#'   probability integral transforms of the data under the margin (iid uniform if
#'   the model is correct), or with `type = "margin"` the probability integral
#'   transforms themselves.
#' }
#'
#' The residuals of [armacop-class], [arfactcop-class] and [sdvinecop-class]
#' models are those of the legacy `tscopula` package, except that tscopula's
#' Gaussian-process residuals are accurate only after a start-up stretch of
#' about ten values, whereas these are exact from the first observation. For a
#' [vtscop-class] model they are not: `tscopula::resid` for a `vtscopula`
#' returns the normal-scale residuals of the base process at \eqn{V(u_t)},
#' which is `resid(model@base, V(u), scale = "normal")` here, whereas the
#' residuals of a [vtscop-class] model are the conditional probability integral
#' transforms of \eqn{u_t} itself, so that they test the model for the
#' observed series rather than for its v-transform.
#'
#' @param object a model with known coefficients or a fitted model object.
#' @param data copula data in (0, 1); for a fitted object it defaults to the
#'   data to which the model was fitted.
#' @param scale `"uniform"` (default) or `"normal"` (the `qnorm` of the
#'   uniform residuals).
#' @param nodes ([vtscop-class] only) number of quadrature nodes for nonlinear
#'   v-transforms.
#' @param ... not used.
#'
#' @return A numeric vector the length of the data.
#' @name resid-method
NULL

#' Prediction for time series copula models
#'
#' Predictive distribution of the copula variable `h` steps ahead given the
#' past data.
#' \itemize{
#'   \item [armacop-class] and [arfactcop-class]: on the normal scale,
#'   \eqn{x_{n+h} = \Phi^{-1}(u_{n+h})} given the past is normal with the mean
#'   and standard deviation returned by `type = "moments"`; the other types
#'   transform this to the uniform scale. A vector of horizons is allowed for
#'   `type = "moments"`.
#'   \item [sdvinecop-class]: one-step-ahead (`h = 1`) distribution given the
#'   last \eqn{K} values, from the forward Rosenblatt functions and their
#'   inverses; the data must be complete. `type` is `"df"`, `"qf"` or
#'   `"dens"`.
#'   \item [vtscop-class]: the base process's conditional law of
#'   \eqn{V_{n+h}} given \eqn{V(u_1), \ldots, V(u_n)} passed through the
#'   stochastic inverse of the v-transform. Horizons `h > 1` are available for
#'   an `arfactcop` base only; a d-vine base has `h = 1`. `type` is `"df"`,
#'   `"qf"` or `"dens"`.
#'   \item [tsmod-class]: the predictive distribution of \eqn{x_{t+h}} on the
#'   original scale given the data \eqn{x_1, \ldots, x_t}. The copula process
#'   gives the predictive distribution of \eqn{u_{t+h}} given the past
#'   probability integral transforms \eqn{F(x_s)}, and the margin maps it to
#'   the original scale: the distribution function is
#'   \eqn{C_{t+h|t}(F(x))}, the quantile function
#'   \eqn{F^{-1}(C^{-1}_{t+h|t}(p))} and the density
#'   \eqn{f(x) c_{t+h|t}(F(x))}. Here `data` are observations on the original
#'   scale and `x` are values of \eqn{x_{t+h}} (probabilities in (0, 1) for
#'   `"qf"`); `type` is `"df"`, `"qf"` or `"dens"`. See [tsrisk()] for the
#'   predictive value-at-risk and expected shortfall.
#' }
#'
#' @param object a model with known coefficients or a fitted model object.
#' @param data past copula data in (0, 1); for a fitted object it defaults to
#'   the data to which the model was fitted.
#' @param h horizon (number of steps ahead).
#' @param type `"moments"` (Gaussian models; default for them) for a data
#'   frame with columns `h`, `mean` and `sd` of the normal score; `"df"` for
#'   the conditional distribution function of \eqn{u_{n+h}} at `x`; `"qf"` for
#'   its quantile function at probabilities `x`; `"dens"` for its density at
#'   `x`. The default is `"df"` for [sdvinecop-class] and [vtscop-class].
#' @param x points in (0, 1) at which to evaluate `"df"`, `"qf"` or `"dens"`
#'   (probabilities for `"qf"`).
#' @param nodes ([vtscop-class] only) number of quadrature nodes for nonlinear
#'   v-transforms.
#' @param ... not used.
#'
#' @return A data frame (`type = "moments"`) or a numeric vector the length of
#'   `x`.
#' @name predict-method
NULL

#' Common methods for time series copula models
#'
#' Methods for extracting coefficients, the covariance matrix and the
#' log-likelihood of time series copula models and fitted model objects, and
#' for printing them. They apply to [armacop-class], [arfactcop-class],
#' [sdvinecop-class], [vtscop-class], [marg-class], [emarg-class], [tsmod-class],
#' [pc-class] and the corresponding fitted classes (`armacopfit`,
#' `arfactcopfit`, `sdvinecopfit`, `vtscopfit`, `margfit`, `tsmodfit`; see
#' [tsfit()]). For [tsmod-class] the coefficients carry the prefixes `margin.`
#' and `copula.`.
#'
#' * `coef()` of a model returns its coefficients: for [armacop-class] the
#'   order `ar`, `ma`, `sar`, `sma` (`NA` for a structure that has not been
#'   fitted); for [arfactcop-class] the rates `lambda1, ...` then the weights
#'   `w1, ...` (the first `m - 1` weights without noise); for
#'   [sdvinecop-class] the coefficients of the kpacf model followed by the
#'   auxiliary parameters of the explicit pair copulas (`nu<lag>` for t
#'   copulas, the second BB1 parameter `delta<lag>` for BB1 copulas) and, with
#'   inverse-v-transformed pair copulas, `delta`, or `delta1` and `delta2`;
#'   for [vtscop-class] those of the base process followed by the free
#'   parameters of the v-transform (`delta`, `kappa`, `xi`). `coef()` of a
#'   fitted object returns the estimates.
#' * `vcov()` of a fitted object returns the covariance matrix of the
#'   estimates (`NA` unless standard errors were requested in `tsfit()`).
#' * `logLik()` of a fitted object returns the maximised log-likelihood with
#'   `df` and `nobs` attributes, so that [stats::AIC()] and [stats::BIC()]
#'   work.
#' For [marg-class] the coefficients are the free parameters (those not in
#' `fix`), in the order of the family's parameters (`mu`, `sigma`, shapes).
#'
#' * `show()` prints a summary of the model or of the fit.
#'
#' @param object a model or fitted model object.
#' @param ... not used.
#'
#' @return `coef()`: a named numeric vector; `vcov()`: a matrix; `logLik()`: an
#'   object of class `"logLik"`; `show()`: the object, invisibly.
#' @name tscop-methods
NULL
