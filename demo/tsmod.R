# Full time series models (tsmod): a margin and a copula process. A series x_t
# with a skew Student t margin and a serial dependence described by an
# AR-factor copula process (and, below, the process-level v-transform model
# for volatile series): simulation, two-stage estimation (the margin by maximum
# likelihood, then the copula process on the probability integral transforms),
# the semiparametric version with an empirical margin, residual diagnostics,
# predictive distributions and risk measures.
# Run with demo("tsmod", package = "tscop") or
#   devtools::load_all(); source("demo/tsmod.R", echo = TRUE)

if (!("package:tscop" %in% search())) library(tscop)
section <- function(title) cat("\n\n==== ", title, " ====\n\n", sep = "")
op <- par(no.readonly = TRUE)
set.seed(1)

section("1. A model and a simulated series")
truth <- tsmod(marg("sst", c(mu = 0, sigma = 1, df = 6, gamma = 1.2)), arfactcop(lambda = 0.9, w = 0.6))
truth
coef(truth)
x <- tssim(truth, 2500)
par(mfrow = c(1, 3))
plot(x, type = "l", main = "simulated series", ylab = "x")
acf(x, main = "ACF")
hist(x, breaks = 50, freq = FALSE, main = "marginal distribution", xlab = "x")
curve(mdens(truth@margin, x), add = TRUE, col = 2, lwd = 2)
par(op)

section("2. Two-stage estimation")
fit <- tsfit(x, tsmod("sst", arfactcop(order = 1)))
fit
rbind(true = coef(truth), estimate = coef(fit))
c(logLik = as.numeric(logLik(fit)), AIC = AIC(fit), BIC = BIC(fit))

section("3. Residual diagnostics")
r <- resid(fit)                         # copula residuals of the PITs: iid uniform if the model is right
par(mfrow = c(1, 3))
hist(r, main = "copula residuals", xlab = "u")
acf(qnorm(r), main = "ACF of normal residuals")
acf(abs(qnorm(r)), main = "ACF of absolute residuals")
par(op)
ks.test(r, "punif")$p.value
plot(x, resid(fit, type = "margin"), pch = 20, cex = 0.4, xlab = "x", ylab = "F(x)",
     main = "fitted margin")

section("3b. One-step-ahead prediction and risk measures")
g <- seq(-4, 5, length.out = 300)
plot(g, predict(fit, type = "dens", x = g), type = "l", xlab = "x", ylab = "density",
     main = "predictive density of the next value")
abline(v = x[length(x)], lty = 3)             # the last observation
predict(fit, type = "qf", x = c(0.05, 0.5, 0.95))   # predictive quantiles
tsrisk(fit, p = c(0.01, 0.05))                # predictive value-at-risk and expected shortfall (lower tail)
tsrisk(fit, p = 0.05, tail = "upper")
tsrisk(fit, p = 0.05, h = 5)                  # five days ahead

section("4. Choosing the margin")
a <- automarg(x)                        # ranks the nine families
a@table[1:4, 1:6]
fa <- tsfit(x, tsmod(a@best@model, arfactcop(order = 1)))   # the selected family as the margin
c(sst = as.numeric(logLik(fit)), selected = as.numeric(logLik(fa)))

section("5. Semiparametric estimation: an empirical margin")
fs <- tsfit(x, tsmod(emarg(x), arfactcop(order = 1)))
fs
round(rbind(parametric = coef(fit)[c("copula.lambda1", "copula.w1")],
            semiparametric = coef(fs)), 3)

section("6. A volatile series: process-level v-transform copula")
vt_true <- tsmod(marg("st", c(mu = 0, sigma = 1, df = 8)),
                 vtscop(arfactcop(lambda = 0.95, w = 0.5), udp::v2p(0.45, 0.8)))
xv <- tssim(vt_true, 2000)
par(mfrow = c(1, 2))
plot(xv, type = "l", main = "series with volatility clustering", ylab = "x")
acf(abs(xv), main = "ACF of absolute values")
par(op)
fv <- tsfit(xv, tsmod("st", vtscop(arfactcop(order = 1), udp::v2p(0.5, 1))))
fv
rbind(true = coef(vt_true), estimate = coef(fv))
ks.test(resid(fv), "punif")$p.value
cat("\nEnd of demo.\n")
