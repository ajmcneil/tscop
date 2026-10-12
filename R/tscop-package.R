#' @references Bladt, M., Dias, A., Han, J. and McNeil, A. J. (2026). Semiparametric and
#'   parametric distributional forecasting of univariate time series using
#'   non-Gaussian ARMA models based on D-vines. *Canadian Journal of Statistics*,
#'   54(1), e70033. \doi{10.1002/cjs.70033}
#'
#'   Dias, A., Han, J. and McNeil, A. J. (2027). GARCH copulas, v-transforms and
#'   D-vines for stochastic volatility. *Journal of Multivariate Analysis*, 217,
#'   105695. \doi{10.1016/j.jmva.2026.105695}
#'
#' @keywords internal
"_PACKAGE"

#' @import methods
#' @importFrom stats coef logLik vcov resid residuals predict pnorm qnorm dnorm optim
#' @importFrom utils modifyList
#' @importFrom rvinecopulib bicop_dist dbicop hbicop ktau_to_par par_to_ktau
#' @importFrom udp pcoincide
NULL
