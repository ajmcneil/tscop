# tscop (development version)

* `show` for a fitted `sdvinecop` now names the base copula of an inverse-v-transformed (`"ivt"`) tail, like `show` for the model.
* `sdvinecop` models with inverse-v-transformed (`"ivt"`) pair copulas accept two-parameter nonlinear v-transforms: `vttype = "v2p"` gives each pair the transforms `v2p(delta1, kappa1)` of the earlier and `v2p(delta2, kappa2)` of the later variable at lags `1, ..., L` (`L` from the kpacf, `nl` and `nltau`), linear ones with the same fulcrums beyond. `vtshape` holds the shape parameters (1 is the linear transform) and `vtfix` the parameters held fixed (any of `delta1`, `delta2`, `kappa1`, `kappa2`; for example `vtfix = "kappa2"`). Fits are staged (linear model first) and store the linear log-likelihood in `loglik_linear`. Density, h-functions, gradient, simulation, residuals and prediction are compiled for the Clayton, Gumbel, Joe, Frank and ast bases (the Gaussian base uses the R implementation). `vtpar` is deprecated: use `vtfix`; `vtpar = "common"` is the process-level model of `vtscop` with a linear v-transform. `autodvinecop` does not accept `vttype = "v2p"`.
* `vtscop` predictive distribution and quantile functions for a nonlinear v-transform are much faster: the integral over the base law is tabulated once per prediction (quantile functions of linear v-transforms are explicit).
* `tsroll()` computes out-of-sample one-step forecasts of a `tsmod`, `sdvinecop`, `vtscop`, `armacop` or `arfactcop` model at every origin of a holdout (predictive quantiles, PIT, log predictive density, lower-tail VaR and ES), with fixed parameters (one compiled pass for `sdvinecop` copulas, including nonlinear v-transforms) or with periodic expanding or rolling re-estimation started from the previous coefficients; the result is a `tsforecast` object. `tsscore()` evaluates one or several sets of forecasts: average quantile score and its profile by level (with a `plot` method), log score, Kupiec and Christoffersen VaR backtests, ES backtests (FZ0, Acerbi-Szekely Z1 and Z2) and Diebold-Mariano tests with Newey-West errors.

# tscop 0.1.0

First version: a rebuild of the `tscopula` package on top of `udp`.

## Copula processes

* Gaussian processes: `armacop` (ARMA and seasonal ARMA) and `arfactcop` (independent AR(1) factors plus noise), with exact Kalman-filter likelihoods.
* `sdvinecop`: stationary d-vine copula processes with a Kendall partial autocorrelation function from an ARMA-type kpacf, any `rvinecopulib` parametric pair copula, BB1, Student t, the `ast` family, and inverse-v-transformed (`"ivt"`) pair copulas for stochastic volatility; automatic selection with `autodvinecop`.
* `vtscop`: process-level v-transforms (any `udp` v-transform) of a Gaussian or d-vine base process.
* `tsmod`: a margin plus a copula process, with two-stage, semiparametric (`emarg`) and joint maximum likelihood estimation, residuals, prediction and risk measures (`tsrisk`).

## Margins

* `marg` families gauss, laplace, slaplace, st, sst, jsu, sas (sinh-arcsinh), nig, hyp and gh; `automarg` chooses among them by AIC, BIC or AICc.

## Inference and diagnostics

* `tsfit` with bootstrap standard errors (multiple cores supported) and exact gradients.
* `AICc` for all fitted classes, alongside `AIC` and `BIC`.
* `ktauhat`: empirical Kendall partial autocorrelations compared with those of a model, with a simulation band.

## Speed

* Compiled (Rcpp) likelihood, gradient, simulation, Rosenblatt transform and prediction kernels for d-vine copulas with all supported pair families, and Kalman kernels for the Gaussian processes. The R implementations remain as references and fallbacks (`options(tscop.cpp = FALSE)`).

## Known limitations

* Rotated `ast` pairs use the slower R path.
* Hessian standard errors are not valid for the fulcrum parameters of IVT pair copulas with Gumbel or Gaussian bases (likelihood cusps); use bootstrap standard errors.
* Nonlinear v-transforms are available at the process level (`vtscop`) but not inside pairwise IVT copulas.
* There are no plotting functions; see the demos for diagnostics.
