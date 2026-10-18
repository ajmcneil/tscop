# Demonstration of inverse-v-transformed (IVT) pair copulas in d-vines, the
# models of Dias, Han and McNeil (2027) for stochastic volatility: each pair
# copula is c(u, v) = c*(V1(u), V2(v); theta_k) with a base copula c* that has
# positive dependence and upper tail dependence (ast, joe, gumbel,
# clayton180) and linear v-transforms V1, V2 with fulcrums delta1, delta2
# shared by all lags; theta_k is calibrated to a kpacf that must be
# non-negative (here an AR-factor kpacf). Run with
#   demo("ivt", package = "tscop")   or
#   devtools::load_all(); source("demo/ivt.R", echo = TRUE)
# It takes about 45 seconds: the likelihood, its exact gradient, the simulation and the
# prediction run in compiled code (see demo("speed", package = "tscop")).

if (!("package:tscop" %in% search())) library(tscop)
set.seed(2029)
ncores <- if (.Platform$OS.type == "windows") 1L else max(1L, parallel::detectCores() - 1L)   # processes for the bootstrap refits
op <- par(no.readonly = TRUE)
section <- function(title) cat("\n\n==== ", title, " ====\n\n", sep = "")

# ---------------------------------------------------------------------------
section("1. An IVT d-vine and its cross-shaped lag-1 dependence")
kp <- arfactcop(lambda = 0.8, w = 0.35)     # ARMA(1,1)-type kpacf, non-negative
mod <- sdvinecop(kp, tail = pc("ivt", base = "ast"), vt = c(0.45, 0.55), tautol = 1e-3)
mod                                         # fulcrums 0.45 (earlier) and 0.55 (later variable)
round(ktau(mod, 6), 3)                      # the kpacf = Kendall's tau of the *base* copulas
u <- tssim(mod, 3000)
gauss <- tssim(sdvinecop(kp, tautol = 1e-3), 3000)
par(mfrow = c(1, 3))
plot(u, type = "l", main = "IVT d-vine", xlab = "t", ylab = "u")
plot(u[-3000], u[-1], pch = 20, cex = 0.4, xlab = "u[t]", ylab = "u[t+1]",
     main = "lag-1: cross shape")
plot(gauss[-3000], gauss[-1], pch = 20, cex = 0.4, xlab = "u[t]", ylab = "u[t+1]",
     main = "Gaussian kpacf-matched")
par(op)

# ---------------------------------------------------------------------------
section("2. The base copulas")
# ast, joe, gumbel and clayton180 have upper tail dependence; frank has none and gauss is the
# cheapest. (Fulcrum estimates are regular for joe, clayton180 and frank; for gumbel and gauss
# the likelihood has cusps at the observations, see ?pc.)
par(mfrow = c(2, 3))
for (b in c("ast", "joe", "gumbel", "clayton180", "frank", "gauss")) {
  m <- sdvinecop(kp, tail = pc("ivt", base = b), vt = c(0.5, 0.5), tautol = 1e-3)
  x <- tssim(m, 2500)
  plot(x[-2500], x[-1], pch = 20, cex = 0.4, xlab = "u[t]", ylab = "u[t+1]", main = b)
}
par(op)

# ---------------------------------------------------------------------------
section("3. Simulate a GARCH(1,1) series and fit the IVT d-vine")
# GARCH(1,1) with t4 innovations (the paper's setting).
sim_garch <- function(n, a0 = 0.1, a1 = 0.3, b1 = 0.6, df = 4, burn = 1000) {
  eps <- rt(n + burn, df) / sqrt(df / (df - 2))
  x <- numeric(n + burn)
  s2 <- a0 / (1 - a1 - b1)
  for (t in seq_along(x)) { x[t] <- sqrt(s2) * eps[t]; s2 <- a0 + a1 * x[t]^2 + b1 * s2 }
  x[-(1:burn)]
}
x <- sim_garch(10000)                         # the length used in the paper
pobs <- rank(x) / (length(x) + 1)             # pseudo-observations
spec <- sdvinecop(arfactcop(order = 1), tail = pc("ivt", base = "ast"), maxlag = 40, tautol = 1e-3)
t0 <- Sys.time()
fit <- tsfit(pobs, spec, pseudo = TRUE)
fit
cat("fitted in", format(round(as.numeric(difftime(Sys.time(), t0, units = "secs")))), "seconds\n")
# The implied ARMA(1,1) kpacf (the paper reports phi = 0.780, psi = -0.614,
# delta1 = 0.510, delta2 = 0.515 for 10000 simulated values):
arfact2arma(fit@model@kpacf)
tau <- ktau(fit, 2)
round(rbind(tau = tau, nu = udp::astcopula_nu(tau)), 3)   # ast degrees of freedom at lags 1 and 2

# ---------------------------------------------------------------------------
section("4. Gaussian alternative and fulcrum options")
# A Gaussian d-vine sees no rank correlation in volatility clustering (large
# moves of either sign follow large moves), so its log-likelihood is about zero
# and its AIC about 4; the IVT d-vine captures the dependence. Fixed (0.5) or
# common fulcrums need fewer parameters when the data are (nearly) symmetric.
gfit <- tsfit(pobs, sdvinecop(arfactcop(order = 1), maxlag = 20, tautol = 1e-3), pseudo = TRUE)
# (one common fulcrum is the process-level model: use vtscop() with a linear v-transform; the deprecated
# vtpar = "common" of sdvinecop() gives the same fit)
cfit <- tsfit(pobs, vtscop(sdvinecop(arfactcop(order = 1), maxlag = 20, tautol = 1e-3, tail = pc("ast")),
                           udp::vlinear(0.5)), pseudo = TRUE)                     # one common fulcrum
sfit <- tsfit(pobs, sdvinecop(arfactcop(order = 1), vtfix = c("delta1", "delta2"), maxlag = 20, tautol = 1e-3,
                             tail = pc("ivt", base = "ast")), pseudo = TRUE)       # fulcrums fixed at 0.5
c(gaussian = AIC(gfit), ivt_both = AIC(fit), ivt_common = AIC(cfit), ivt_fixed = AIC(sfit))

# ---------------------------------------------------------------------------
section("4b. Empirical Kendall partial autocorrelations")
# For IVT models the generalised-lagged pairs are v-transformed, so the diagnostic is on the scale
# of the base copulas, where the kpacf lives; the band is simulated from the fitted model.
kh <- ktauhat(fit, lagmax = 30, B = 50, cores = ncores)
round(kh[c(1:5, 10, 20, 30), ], 3)

# ---------------------------------------------------------------------------
section("5. Diagnostics and one-step predictive distributions")
r <- resid(fit, scale = "normal")
par(mfrow = c(1, 3))
qqnorm(r, main = "residuals"); qqline(r)
acf(r, main = "ACF of residuals")
acf(abs(r), main = "ACF of absolute residuals")
par(op)
# After a calm day and after a large move the predictive density of the next
# day's return differs sharply: large moves of either sign are likely after a
# large move (the cross shape), but the model also reproduces the asymmetry
# through the fulcrums.
grid <- seq(0.005, 0.995, length.out = 199)
calm <- c(pobs[-length(pobs)], 0.5)
shock <- c(pobs[-length(pobs)], 0.03)
plot(grid, predict(fit, calm, type = "dens", x = grid), type = "l", lwd = 2, ylim = c(0, 3),
     xlab = "next value (copula scale)", ylab = "predictive density",
     main = "one-step predictive density")
lines(grid, predict(fit, shock, type = "dens", x = grid), lwd = 2, col = 2)
legend("top", c("after a median value", "after a value in the lower 3%"), lwd = 2, col = 1:2, bty = "n")
# One-step quantiles on the copula scale (they would be mapped through a margin
# to give value-at-risk):
predict(fit, type = "qf", x = c(0.01, 0.05, 0.5, 0.95, 0.99))

# ---------------------------------------------------------------------------
section("6. Standard errors")
fit_h <- tsfit(pobs, spec, pseudo = TRUE, se = "hessian")
fit_h
# The parametric bootstrap (the recommended standard errors) with 100 replicates: each one is
# simulated from the fitted model, re-ranked (pseudo = TRUE) and refitted from the fitted values.
system.time(fit_b <- tsfit(pobs, spec, pseudo = TRUE, se = "bootstrap", B = 100, cores = ncores))
rbind(hessian = fit_h@se, bootstrap = fit_b@se)
cat("\nEnd of demo.\n")
