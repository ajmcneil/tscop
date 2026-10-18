# Out-of-sample forecasting and forecast evaluation: tsroll() computes, at every origin of a holdout, the
# one-step predictive distribution of the next observation (quantiles, probability integral transform,
# log density, lower-tail VaR and ES), with fixed parameters or periodic re-estimation; tsscore()
# scores and backtests one or several sets of forecasts. Run with
#   demo("forecast", package = "tscop")   or
#   devtools::load_all(); source("demo/forecast.R", echo = TRUE)
# It takes about 2-3 minutes (the fits dominate; the holdout passes run in compiled code).

if (!("package:tscop" %in% search())) library(tscop)
set.seed(2033)
op <- par(no.readonly = TRUE)
section <- function(title) cat("\n\n==== ", title, " ====\n\n", sep = "")

# ---------------------------------------------------------------------------
section("1. A volatile series with nonlinear IVT dependence")
truth <- tsmod(marg("st", c(mu = 0, sigma = 0.01, df = 5)),
               sdvinecop(arfactcop(lambda = 0.9, w = 0.3), tail = pc("ivt", base = "joe"),
                         vt = c(0.5, 0.55), vttype = "v2p", vtshape = c(2, 0.6), maxlag = 20, tautol = 1e-3))
x <- tssim(truth, 1500)
plot(x, type = "l", xlab = "t", ylab = "return", main = "simulated series")
abline(v = 1000, lty = 3)           # estimation sample: 1-1000, holdout: 1001-1500

# ---------------------------------------------------------------------------
section("2. Fit competing models on the first 1000 observations")
spec <- function(...) tsmod(marg("st"), sdvinecop(arfactcop(order = 1), tail = pc("ivt", base = "joe"),
                                                  maxlag = 20, tautol = 1e-3, ...))
fit_lin <- tsfit(x[1:1000], spec())
fit_nl <- tsfit(x[1:1000], spec(vttype = "v2p"))
c(linear = AIC(fit_lin), v2p = AIC(fit_nl))

# ---------------------------------------------------------------------------
section("3. Fixed-parameter holdout forecasts")
# With refit = Inf (the default) the parameters stay at the fitted values while each forecast conditions
# on all the data before its origin. The true model is the benchmark that cannot be beaten systematically.
fc_lin <- tsroll(fit_lin, x, start = 1001)
fc_nl <- tsroll(fit_nl, x, start = 1001)
fc_true <- tsroll(truth, x, start = 1001)
fc_nl                                                # forecasts and VaR breaches at p = 1% and 5%
sc <- tsscore(fc_lin, v2p = fc_nl, true = fc_true)
sc                                                   # AQS, log score, VaR and ES backtests, Diebold-Mariano tests
# The AQS profile shows at which levels the models differ; negative values favour the model over the first.
plot(sc, relative = TRUE, main = "quantile scores relative to the linear model")

# ---------------------------------------------------------------------------
section("4. Periodic re-estimation")
# Every 250 observations the model is re-estimated on an expanding window, starting from the previous
# coefficients; refit = H with window = "rolling" and width = w uses the last w observations instead.
fc_re <- tsroll(fit_nl, x, start = 1001, refit = 250)
fc_re@blocks                                         # window and status of each block
round(fc_re@estimates, 3)                            # coefficients of the model of each block
tsscore(fc_nl, refit = fc_re)

# ---------------------------------------------------------------------------
section("5. Forecast risk measures and calibration")
t <- fc_nl@index
plot(t, fc_nl@observed, type = "l", col = "grey50", xlab = "t", ylab = "return", main = "returns and 5% VaR / ES")
lines(t, fc_nl@var[, 2], col = 2, lwd = 2)
lines(t, fc_nl@es[, 2], col = 4, lwd = 2)
breach <- fc_nl@observed < fc_nl@var[, 2]
points(t[breach], fc_nl@observed[breach], pch = 20, col = 2)
legend("bottomleft", c("VaR 5%", "ES 5%", "breach"), col = c(2, 4, 2), lwd = c(2, 2, NA), pch = c(NA, NA, 20), bty = "n")
# The PIT values of a correct model are iid uniform:
par(mfrow = c(1, 2))
hist(fc_nl@pit, breaks = 20, freq = FALSE, main = "PIT", xlab = "u"); abline(h = 1, lty = 2)
acf(abs(qnorm(fc_nl@pit)), main = "ACF of |normal scores|")
par(op)
# Forecasts from other sources can be scored too, as matrices of predictive quantiles:
tsscore(fc_nl@quantile, observed = fc_nl@observed, levels = fc_nl@levels, p = 0.05)
cat("\nEnd of demo.\n")
