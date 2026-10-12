# Process-level v-transform copula models (vtscop): a base copula process of
# non-negative Kendall partial autocorrelations, V_t, is turned into a copula
# process for U_t by the stochastic inverse of a v-transform. The likelihood is
# the base likelihood at V(u_t), for ANY v-transform, so the nonlinear
# families of the udp package (v2p, v3p, v3b, ...) are fitted in seconds.
# Takes about half a minute. Run with demo("vtscop", package = "tscop") or
#   devtools::load_all(); source("demo/vtscop.R", echo = TRUE)

if (!("package:tscop" %in% search())) library(tscop)
section <- function(title) cat("\n\n==== ", title, " ====\n\n", sep = "")
op <- par(no.readonly = TRUE)
set.seed(1)

section("1. Specification and simulation")
# One AR factor plus noise as base process, asymmetric 2-piece v-transform
m <- vtscop(arfactcop(lambda = 0.95, w = 0.5), udp::v2p(delta = 0.45, kappa = 0.8))
m
coef(m)
u <- tssim(m, 3000)
par(mfrow = c(1, 3))
plot(u[1:500], type = "l", main = "simulated copula series", ylab = "u")
acf(abs(qnorm(u)), main = "ACF of |qnorm(u)| (volatility clustering)")
plot(u[-length(u)], u[-1], pch = 20, cex = 0.3, xlab = "u[t]", ylab = "u[t+1]",
     main = "lag-1 scatter")
par(op)

section("2. Fitting, with standard errors")
fit <- tsfit(u, vtscop(arfactcop(order = 1), udp::v2p(0.5, 1)), pseudo = FALSE, se = "hessian")
fit
rbind(true = coef(m), estimate = coef(fit))
# Fixing the fulcrum at its true value
fit2 <- tsfit(u, vtscop(arfactcop(order = 1), udp::v2p(0.45, 1), vtfix = "fulcrum"), pseudo = FALSE)
coef(fit2)

section("3. Residuals and prediction")
r <- resid(fit, scale = "normal")
par(mfrow = c(1, 3))
qqnorm(r, main = "residuals"); qqline(r)
acf(r, main = "ACF of residuals")
acf(abs(r), main = "ACF of absolute residuals")
par(op)
# Predictive density of the next value, one and five days ahead
x <- seq(0.005, 0.995, length.out = 400)
plot(x, predict(fit, u, type = "dens", x = x), type = "l", ylab = "predictive density",
     xlab = "u", main = "predictive density of the next copula value")
lines(x, predict(fit, u, h = 5, type = "dens", x = x), lty = 2)
legend("top", c("h = 1", "h = 5"), lty = 1:2, bty = "n")
predict(fit, u, type = "qf", x = c(0.05, 0.5, 0.95))

cat("\nEnd of demo.\n")

