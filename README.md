# tscop

A ground-up rebuild of the `tscopula` package for time series copula models, built on [udp](https://github.com/ajmcneil/udp). The legacy `tscopula` package is left untouched.

## What is in the package

- **Copula processes**: `armacop` and `arfactcop` (Gaussian), `sdvinecop` (stationary d-vine processes with an ARMA-type Kendall partial autocorrelation function, including non-Gaussian pair copulas chosen automatically by `autodvinecop`, and inverse-v-transformed pair copulas, with linear or nonlinear v-transforms, for stochastic volatility), and `vtscop` (process-level v-transforms).
- **Margins**: `marg` (Gaussian, Laplace, Student t, skew t, Johnson SU, sinh-arcsinh, NIG, hyperbolic, generalised hyperbolic), `automarg`, and `emarg` for the semiparametric approach.
- **Models and inference**: `tsmod` combines a margin and a copula process; `tsfit` fits by two-stage or joint maximum likelihood with Hessian or bootstrap standard errors; `AIC`, `BIC` and `AICc` compare models; `tssim`, `ktauhat`, `tsrisk` and the `plot` method for pair copulas simulate, diagnose and display the dependence structure.
- **Forecasting**: `tsroll` computes one-step-ahead predictive distributions over a holdout, with fixed or re-estimated parameters; `tsscore` evaluates them (quantile score, log score, risk-measure backtests, Diebold-Mariano tests) and `tstable` tabulates the results.
- **Speed**: compiled likelihood, exact gradient, simulation and prediction kernels (Rcpp) for all d-vine pair copula families.
- **Data**: `indpro`, monthly US industrial production (Federal Reserve Board).

## Getting started

Two vignettes show the main workflow from data to forecast evaluation:

- `vignette("ip", package = "tscop")`: non-Gaussian ARMA dependence for industrial production growth, with `sdvinecop` models with a Gaussian tail, `autodvinecop`, `tsroll`, `tsscore` and `tstable`.
- `vignette("vix", package = "tscop")`: nonlinear inverse-v-transformed d-vine models for VIX returns, with semiparametric, IFM and joint estimation, compared with GARCH models.

Shorter worked examples are in `demo()`: `demo("armacop", package = "tscop")`, and likewise `arfactcop`, `sdvinecop`, `ivt`, `vtscop`, `marg`, `tsmod`, `speed`, `forecast` and `ivt-nonlinear`. See `NEWS.md` for the changes in each version.

## Installation

```r
# install.packages("remotes")
remotes::install_github("ajmcneil/tscop", build_vignettes = TRUE)
```

The [udp](https://github.com/ajmcneil/udp) package, which is installed from GitHub, is pulled in automatically through the `Remotes` field of `DESCRIPTION`.

## License

GPL-3; see `LICENSE`.

## Background and references

The package implements and extends models from the following papers (BibTeX entries in `inst/REFERENCES.bib`):

- McNeil, A. J. (2021). Modelling volatile time series with v-transforms and copulas. *Risks*, 9(1), 14. [doi:10.3390/risks9010014](https://doi.org/10.3390/risks9010014). V-transforms for volatile financial series, with a Gaussian ARMA copula process as the base process.
- Bladt, M. and McNeil, A. J. (2022). Time series copula models using d-vines and v-transforms. *Econometrics and Statistics*, 24, 27-48. [doi:10.1016/j.ecosta.2021.07.004](https://doi.org/10.1016/j.ecosta.2021.07.004). V-transforms combined with non-Gaussian d-vine base processes.
- Bladt, M. and McNeil, A. J. (2022). Time series with infinite-order partial copula dependence. *Dependence Modeling*, 10(1), 87-107. [doi:10.1515/demo-2022-0105](https://doi.org/10.1515/demo-2022-0105). Models parameterised through the Kendall partial autocorrelation function.
- Bladt, M., Dias, A., Han, J. and McNeil, A. J. (2026). Semiparametric and parametric distributional forecasting of univariate time series using non-Gaussian ARMA models based on D-vines. *Canadian Journal of Statistics*, 54(1), e70033. [doi:10.1002/cjs.70033](https://doi.org/10.1002/cjs.70033). Non-Gaussian ARMA models based on d-vines (`sdvinecop`).
- Dias, A., Han, J. and McNeil, A. J. (2027). GARCH copulas, v-transforms and D-vines for stochastic volatility. *Journal of Multivariate Analysis*, 217, 105695. [doi:10.1016/j.jmva.2026.105695](https://doi.org/10.1016/j.jmva.2026.105695). Inverse-v-transformed pair copulas and the `ast` family (the copula of the absolute values of a spherical t distribution).
