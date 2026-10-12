# tscop

A ground-up rebuild of the `tscopula` package for time series copula models, built on [udp](https://github.com/ajmcneil/udp). Under development; the legacy `tscopula` package is left untouched.

## Installation

```r
# install.packages("remotes")
remotes::install_github("ajmcneil/tscop")
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
