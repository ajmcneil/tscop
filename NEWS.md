# tscop (development version)

* `show` for a fitted `sdvinecop` now names the base copula of an inverse-v-transformed (`"ivt"`) tail, like `show` for the model.
* `vtscop` predictive distribution and quantile functions for a nonlinear v-transform are much faster: the integral over the base law is tabulated once per prediction (quantile functions of linear v-transforms are explicit).

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
