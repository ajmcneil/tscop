# Nonlinear (v2p) v-transforms in the inverse-v-transformed (IVT) d-vines of demo("ivt"). Each pair
# copula at lags 1, ..., L is c*(T1(u), T2(v); theta_k), where the v-transforms T1, T2 are the
# two-parameter udp::v2p transforms with fulcrums delta1, delta2 and shapes kappa1, kappa2 (kappa = 1
# is the linear transform of demo("ivt")). Beyond lag L, which is read off the kpacf, the transforms
# are linear with the same fulcrums. Run with
#   demo("ivt-nonlinear", package = "tscop")   or
#   devtools::load_all(); source("demo/ivt-nonlinear.R", echo = TRUE)
# It takes about 1-2 minutes: the likelihood, its exact gradient, the simulation and the prediction
# run in compiled code.

if (!("package:tscop" %in% search())) library(tscop)
set.seed(2031)
op <- par(no.readonly = TRUE)
section <- function(title) cat("\n\n==== ", title, " ====\n\n", sep = "")

# ---------------------------------------------------------------------------
section("1. A v2p IVT d-vine")
kp <- arfactcop(lambda = 0.9, w = 0.3)            # non-negative kpacf
lin <- sdvinecop(kp, tail = pc("ivt", base = "joe"), vt = c(0.5, 0.55), tautol = 1e-3)
nlm <- sdvinecop(kp, tail = pc("ivt", base = "joe"), vt = c(0.5, 0.55), tautol = 1e-3,
                 vttype = "v2p", vtshape = c(2, 0.6))   # kappa1 = 2 (earlier), kappa2 = 0.6 (later variable)
nlm                                              # lags with nonlinear transforms: L, from the kpacf
# The transforms at the fulcrums (0.5, 0.55): kappa bends the branches of the V shape.
v <- seq(0, 1, length.out = 201)
par(mfrow = c(1, 2))
plot(v, udp::udptrans(udp::v2p(0.5, 2), v), type = "l", lwd = 2, xlab = "u", ylab = "T1(u)", main = "earlier variable (kappa1 = 2)")
abline(0, 1, lty = 3)
plot(v, udp::udptrans(udp::v2p(0.55, 0.6), v), type = "l", lwd = 2, xlab = "v", ylab = "T2(v)", main = "later variable (kappa2 = 0.6)")
abline(0, 1, lty = 3)
par(op)

# ---------------------------------------------------------------------------
section("2. Simulation: the shape of the lag-1 dependence")
n <- 3000
x_lin <- tssim(lin, n)
x_nl <- tssim(nlm, n)
par(mfrow = c(1, 2))
plot(x_lin[-n], x_lin[-1], pch = 20, cex = 0.4, xlab = "u[t]", ylab = "u[t+1]", main = "linear (kappa = 1)")
plot(x_nl[-n], x_nl[-1], pch = 20, cex = 0.4, xlab = "u[t]", ylab = "u[t+1]", main = "v2p (kappa = 2, 0.6)")
par(op)
# The kpacf (Kendall's tau of the base copulas) is the same in both models:
round(rbind(linear = ktau(lin, 5), v2p = ktau(nlm, 5)), 3)

# ---------------------------------------------------------------------------
section("3. Fitting: staged estimation")
# The v2p fit starts from the linear model (fulcrums only), then releases the shapes; the linear
# log-likelihood is stored. The number of lags with nonlinear transforms is held fixed within an
# optimisation and updated between refits.
x <- tssim(nlm, 1500)
spec_lin <- sdvinecop(arfactcop(order = 1), tail = pc("ivt", base = "joe"), maxlag = 20, tautol = 1e-3)
spec_nl <- sdvinecop(arfactcop(order = 1), tail = pc("ivt", base = "joe"), maxlag = 20, tautol = 1e-3,
                     vttype = "v2p")
t0 <- Sys.time()
fit_lin <- tsfit(x, spec_lin, pseudo = TRUE)
fit_nl <- tsfit(x, spec_nl, pseudo = TRUE)
cat("fitted in", format(round(as.numeric(difftime(Sys.time(), t0, units = "secs")))), "seconds\n")
fit_nl
c(linear = as.numeric(logLik(fit_lin)), v2p = as.numeric(logLik(fit_nl)),
  v2p_stage1 = fit_nl@loglik_linear)              # stage 1 = the linear fit with the same structure
c(linear = AIC(fit_lin), v2p = AIC(fit_nl))       # the shapes cost two parameters
# The data were simulated with fulcrums (0.5, 0.55) and shapes (2, 0.6).

# ---------------------------------------------------------------------------
section("4. Holding parameters fixed")
# vtfix fixes any of delta1, delta2, kappa1, kappa2 at their values in the specification; for example
# a nonlinear transform for the later variable only:
spec_k2 <- sdvinecop(arfactcop(order = 1), tail = pc("ivt", base = "joe"), maxlag = 20, tautol = 1e-3,
                     vttype = "v2p", vtshape = c(1, 1), vtfix = "kappa1")
fit_k2 <- tsfit(x, spec_k2, pseudo = TRUE)
c(both = AIC(fit_nl), kappa2_only = AIC(fit_k2), linear = AIC(fit_lin))

# ---------------------------------------------------------------------------
section("5. Diagnostics and prediction")
r <- resid(fit_nl, scale = "normal")
par(mfrow = c(1, 3))
qqnorm(r, main = "residuals"); qqline(r)
acf(r, main = "ACF of residuals")
acf(abs(r), main = "ACF of absolute residuals")
par(op)
grid <- seq(0.005, 0.995, length.out = 199)
calm <- c(x[-length(x)], 0.5)
shock <- c(x[-length(x)], 0.03)
plot(grid, predict(fit_nl, calm, type = "dens", x = grid), type = "l", lwd = 2, ylim = c(0, 3),
     xlab = "next value (copula scale)", ylab = "predictive density", main = "one-step predictive density")
lines(grid, predict(fit_nl, shock, type = "dens", x = grid), lwd = 2, col = 2)
lines(grid, predict(fit_lin, shock, type = "dens", x = grid), lwd = 2, col = 2, lty = 2)
legend("top", c("v2p, after a median value", "v2p, after a value in the lower 3%", "linear, after a value in the lower 3%"),
       lwd = 2, col = c(1, 2, 2), lty = c(1, 1, 2), bty = "n")
cat("\nEnd of demo.\n")
