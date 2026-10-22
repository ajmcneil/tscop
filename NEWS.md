# tscop (development version)

* `tsscore()` no longer stops with an error when a log predictive density is infinite (for example the kernel density of an empirical margin at an observation far outside the data): the Diebold-Mariano test of the negative log score is then `NA`.
* `tsscore()` with `p` a subset of the probabilities of the forecasts (for example `p = 0.05` for forecasts with `p = c(0.01, 0.05)`) used the value-at-risk and expected shortfall of the first columns instead of the matching ones, so the backtests were labelled with the wrong probability. It now selects the matching columns and stops with an error if the forecasts do not contain a requested probability.
* New vignette `vignette("vix")`: modelling VIX returns with semiparametric, IFM and joint-ML estimation of nonlinear IVT d-vine copula models, compared in sample and out of sample with GARCH models. The vignette is precomputed (`vignettes/vix.Rmd.orig`, knitted by `vignettes/knit-vignettes.R`); it downloads the prices with `quantmod` when the source is knitted and no data are distributed with the package.
* `show` describes an AR-factor copula process (as the kpacf of a d-vine, or on its own) as `ARFACT(m)` with `m` the number of factors, parallel to `ARMA(p,q)` for `armacop`, instead of "m-factor AR(1) model".
* Two-stage fits of nonlinear (v2p) `sdvinecop` models no longer stall right after the staged linear fit. When a kpacf parameter of the linear fit sits on the edge of its feasible range (for example a persistence close to 1), every step of BFGS along the negative gradient is infeasible, so the optimiser stopped at the linear fit (`kappa1 = kappa2 = 1` exactly, the gradient far from zero) while reporting convergence; observed on a real volatility series (order 2, ast base: log-likelihood 147.9 instead of 159.4). If the exact gradient at a BFGS result is not close to zero, the fit now continues with box-constrained L-BFGS-B and, if needed, the remedy of the joint fit (coordinate probes, restarts on a coarser scale, Nelder-Mead). The `convergence` slot is 10 if the fit is still stalled after that. Fits that converged are not slowed (one gradient evaluation per optimisation is added).
* `plot()` for `sdvinecop` models, fitted `sdvinecop` models and fitted `tsmod` models with a d-vine copula draws the pair copulas at chosen lags (`lags`, default 1:4) as filled contour plots in the style of udp's `bsicopula` plot: the partial copula of `(U_t, U_{t+k})` given the observations in between. `data` overlays the data pairs of each lag (conditioned through the h-functions of the lower lags), `strips = TRUE` shows the two v-transforms of inverse-v-transformed pair copulas in marginal strips with the fulcrums marked, `fill = FALSE` gives contour lines only, and `plot(x, y)` with two models draws them side by side (one row per lag) with common contour levels.
* Joint maximum likelihood (`tsfit(method = "joint")`) is robust against stalling and never ends below its start. (1) The start is the two-stage fit, or the parameters of the copula model supplied in `object` (with the margin refitted) if those have the higher likelihood: refitting a nonlinear (v2p) copula from scratch can end in a worse local optimum, and the joint fit then stalled there (observed on real return series: 6 of 144 joint fits ended 1.2 to 11.3 log-likelihood units below the supplied two-stage fit). (2) The optimiser probes the BFGS result with coordinate steps of several sizes, because the objective is rough in the margin location when observations cross the fulcrums of an inverse v-transform; if a probe improves the log-likelihood by 0.01 or more, BFGS is restarted on a coarser scale, with Nelder-Mead as the last resort. The result is never worse than the start. (3) The number of nonlinear lags `L` follows the two-stage rule (a refit only if the log-likelihood changes by 0.01 or more with the new `L`, kept only if it improves the likelihood). (4) `tsmodfit` has new slots `convergence` (0 converged, 1 iteration limit, 10 stalled) and `message`; a fit that did not converge warns and `show` reports it.
* Joint maximum likelihood (`tsfit(method = "joint")`) treats margin parameters that put a PIT within 1e-9 of 0 or 1, or the median PIT outside (0.05, 0.95), as infeasible: the copula likelihood is unbounded when all the PITs collapse onto the boundary, and the optimiser could escape there (observed with a Johnson SU margin and a nonlinear IVT copula).
* The Gaussian base of nonlinear (v2p) `sdvinecop` models is compiled like the other bases (closed-form h-functions and density on the normal-score scale): a cold v2p fit with a Gaussian IVT base at n = 1500 takes about 2 s instead of about 58 s. Results agree with the R reference to 1e-8.
* Faster nonlinear (v2p) `sdvinecop` fits: the optimiser uses relative tolerance `1e-8` (instead of `1e-10`) and parameters scaled to the size of their standard errors (`control$pscale`, default 6 divided by `sqrt(n)`), which roughly halves the number of likelihood evaluations of a cold fit (Joe, n = 1500: about 14 s to 6-8 s) and makes a warm refit (`init` = fitted model) almost free. Estimates change by about 1e-4.
* `tsfit(se = "hessian")` for `sdvinecop` models uses central differences of the exact gradient instead of the value, which is about four times faster for nonlinear (v2p) models.
* The fitting rule for the number of nonlinear lags `L` (v2p) refits only if the log-likelihood at the optimum changes by 0.01 or more with the new `L`, like the rule for `K`.
* `show` for `sdvinecop` models and fits with `vttype = "v2p"` gives the number of nonlinear lags `L` and what determined it (`nltau`, `nl`) on its own line under the effective lag.
* `tsscore()` reports the negative log score (`neglogscore`) so that every score in the summary is a loss (smaller is better); the printed output says so. Diebold-Mariano comparisons are unchanged.
* `show` for a fitted `sdvinecop` now names the base copula of an inverse-v-transformed (`"ivt"`) tail, like `show` for the model.
* `sdvinecop` models with inverse-v-transformed (`"ivt"`) pair copulas accept two-parameter nonlinear v-transforms: `vttype = "v2p"` gives each pair the transforms `v2p(delta1, kappa1)` of the earlier and `v2p(delta2, kappa2)` of the later variable at lags `1, ..., L` (`L` from the kpacf, `nl` and `nltau`), linear ones with the same fulcrums beyond. `vtshape` holds the shape parameters (1 is the linear transform) and `vtfix` the parameters held fixed (any of `delta1`, `delta2`, `kappa1`, `kappa2`; for example `vtfix = "kappa2"`). Fits are staged (linear model first) and store the linear log-likelihood in `loglik_linear`. Density, h-functions, gradient, simulation, residuals and prediction are compiled for the Gaussian, Clayton, Gumbel, Joe, Frank and ast bases. `vtpar` is deprecated: use `vtfix`; `vtpar = "common"` is the process-level model of `vtscop` with a linear v-transform. `autodvinecop` does not accept `vttype = "v2p"`.
* New demos `demo("ivt-nonlinear")` (v2p v-transforms: simulation, staged fitting, fixed parameters, prediction) and `demo("forecast")` (`tsroll()` and `tsscore()`). `tsscore()` labels unnamed forecasts with their argument names.
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
