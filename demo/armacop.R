# Demonstration of the armacop class: Gaussian ARMA / seasonal ARMA copula
# processes. Run it with
#   demo("armacop", package = "tscop")      (after installing tscop), or
#   devtools::load_all(); source("demo/armacop.R", echo = TRUE)
# from the source folder. It produces a few plots.

if (!("package:tscop" %in% search())) library(tscop)
set.seed(2026)
op <- par(no.readonly = TRUE)

section <- function(title) cat("\n\n==== ", title, " ====\n\n", sep = "")

# ---------------------------------------------------------------------------
section("1. Specifying models")
# A known model: ARMA(1,1) with ar = 0.7 and ma = 0.4 (stats::arima sign
# convention). The innovation variance is not a parameter: the normal scores
# qnorm(u_t) have mean zero and variance one, so only the coefficients matter.
mod <- armacop(ar = 0.7, ma = 0.4)
mod
coef(mod)

# Seasonal terms take a period (12 for monthly, 4 for quarterly data).
smod <- armacop(ar = 0.4, ma = 0.3, sar = 0.5, sma = -0.3, period = 12)
smod

# A structure to be fitted: orders but no coefficients.
armacop(order = c(1, 1))
armacop(order = c(1, 0), seasonal = c(1, 1), period = 4)

# Inadmissible models are rejected.
try(armacop(ar = 1.2))
try(armacop(ma = -1.5))

# ---------------------------------------------------------------------------
section("2. Kendall partial autocorrelation function (ktau)")
# ktau() gives Kendall's tau of the pair copulas in the d-vine representation
# at lags 1, 2, ...: (2/pi) asin of the partial autocorrelations.
round(ktau(mod, 8), 3)
round(ktau(smod, 14), 3)
par(mfrow = c(1, 2))
barplot(ktau(mod, 15), names.arg = 1:15, xlab = "lag", ylab = "Kendall's tau",
        main = "kpacf of ARMA(1,1)")
barplot(ktau(smod, 30), names.arg = 1:30, xlab = "lag", ylab = "Kendall's tau",
        main = "kpacf of seasonal model")
par(op)

# ---------------------------------------------------------------------------
section("3. Simulation (tssim)")
# Exact stationary simulation: no burn-in is needed.
u <- tssim(mod, 500)
head(u)
par(mfrow = c(1, 3))
plot(u, type = "l", xlab = "t", ylab = "u", main = "simulated copula data")
plot(u[-500], u[-1], xlab = "u[t]", ylab = "u[t+1]", main = "lag-1 scatter", pch = 20, cex = 0.5)
acf(qnorm(u), main = "ACF of normal scores")
par(op)

# ---------------------------------------------------------------------------
section("4. Fitting (tsfit) to simulated data")
# Data are copula data in (0, 1): pseudo-observations (ranks) or probability
# integral transforms. Here we start from the structure only: starting values
# are computed automatically (Hannan-Rissanen).
u <- tssim(mod, 1000)
fit <- tsfit(u, armacop(order = c(1, 1)), pseudo = FALSE)
fit
logLik(fit)
AIC(fit)
BIC(fit)

# Compare orders by AIC
orders <- list(c(1, 0), c(2, 0), c(1, 1), c(2, 1), c(0, 2))
aic <- sapply(orders, function(o) AIC(tsfit(u, armacop(order = o), pseudo = FALSE)))
names(aic) <- sapply(orders, function(o) paste0("ARMA(", o[1], ",", o[2], ")"))
round(aic, 1)

# ---------------------------------------------------------------------------
section("5. Standard errors")
# se = "hessian": inverse observed information carried to the coefficients by
# the delta method (quick but only approximate for pseudo-observations).
fit_h <- tsfit(u, armacop(order = c(1, 1)), pseudo = FALSE, se = "hessian")
fit_h
vcov(fit_h)

# se = "bootstrap": parametric bootstrap, re-ranking the simulated data when
# pseudo = TRUE (B = 200 replicates, which takes a second or two).
system.time(fit_b <- tsfit(u, armacop(order = c(1, 1)), pseudo = FALSE, se = "bootstrap", B = 200))
fit_b
rbind(hessian = fit_h@se, bootstrap = fit_b@se)

# ---------------------------------------------------------------------------
section("6. Pseudo-observations from non-normal data")
# Real data have their own margin; the copula is fitted to the ranks.
x <- qt(u, df = 5) * 2 + 10 # t(5) margin, location 10, scale 2
pobs <- rank(x) / (length(x) + 1)
fit_p <- tsfit(pobs, armacop(order = c(1, 1)), pseudo = TRUE, se = "hessian")
fit_p

# ---------------------------------------------------------------------------
section("7. Residual diagnostics (resid)")
# Standardised one-step prediction errors: iid U(0,1) (or N(0,1)) under the model.
r <- resid(fit, scale = "normal")
par(mfrow = c(1, 3))
acf(r, main = "ACF of residuals")
qqnorm(r, main = "normal QQ of residuals"); qqline(r)
hist(resid(fit), main = "uniform residuals", xlab = "", freq = FALSE); abline(h = 1, col = 2)
par(op)
# A misspecified model leaves structure in the residuals
bad <- tsfit(u, armacop(order = c(1, 0)), pseudo = FALSE)
c(ARMA11 = AIC(fit), AR1 = AIC(bad))
Box.test(resid(fit, scale = "normal"), lag = 10, type = "Ljung-Box")$p.value
Box.test(resid(bad, scale = "normal"), lag = 10, type = "Ljung-Box")$p.value

# ---------------------------------------------------------------------------
section("8. Prediction (predict)")
# Conditional mean and sd of the normal score h steps ahead:
predict(fit, h = 1:5)
# The one-step predictive distribution of the copula variable:
x <- seq(0.01, 0.99, length.out = 99)
par(mfrow = c(1, 3))
plot(x, predict(fit, type = "dens", x = x), type = "l", xlab = "u", ylab = "density",
     main = "predictive density"); abline(h = 1, lty = 2)
plot(x, predict(fit, type = "df", x = x), type = "l", xlab = "u", ylab = "df",
     main = "predictive df"); abline(0, 1, lty = 2)
plot(x, predict(fit, type = "qf", x = x), type = "l", xlab = "p", ylab = "quantile",
     main = "predictive quantile function"); abline(0, 1, lty = 2)
par(op)
# Predictive 5%, 50% and 95% quantiles for the next observation:
predict(fit, type = "qf", x = c(0.05, 0.5, 0.95))

# ---------------------------------------------------------------------------
section("9. A seasonal (monthly) model")
smod <- armacop(ar = 0.4, sar = 0.6, sma = -0.3, period = 12)
us <- ts(tssim(smod, 360), start = c(1995, 1), frequency = 12)
sfit <- tsfit(us, armacop(order = c(1, 0), seasonal = c(1, 1), period = 12), pseudo = FALSE,
              se = "hessian")
sfit
# The time series index is carried through to the residuals
head(resid(sfit))
plot(resid(sfit, scale = "normal"), main = "residuals of the seasonal fit", ylab = "")

# ---------------------------------------------------------------------------
section("10. Missing values and warm starts")
umiss <- u
umiss[c(10, 11, 200, 500)] <- NA
tsfit(umiss, armacop(order = c(1, 1)), pseudo = FALSE)@nobs
# Starting from an earlier fit (fixed structure, new data) is cheap:
tsfit(tssim(mod, 1000), fit@model, pseudo = FALSE)
# or supply starting values explicitly
tsfit(u, armacop(order = c(1, 1)), init = c(0.3, 0.1), pseudo = FALSE)

# ---------------------------------------------------------------------------
section("11. Speed")
big <- tssim(smod, 3000)
system.time(tsfit(big, armacop(order = c(1, 0), seasonal = c(1, 1), period = 12), pseudo = FALSE))
cat("\nEnd of demo.\n")
