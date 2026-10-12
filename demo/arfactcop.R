# Demonstration of the arfactcop class: Gaussian copula processes made of
# independent AR(1) "factors" plus white noise, whose partial
# autocorrelations are non-negative. Run it with
#   demo("arfactcop", package = "tscop")      (after installing tscop), or
#   devtools::load_all(); source("demo/arfactcop.R", echo = TRUE)
# from the source folder. The theory is in docs/factor-kpacf-theory.Rmd.

if (!("package:tscop" %in% search())) library(tscop)
set.seed(2028)
op <- par(no.readonly = TRUE)
section <- function(title) cat("\n\n==== ", title, " ====\n\n", sep = "")

# ---------------------------------------------------------------------------
section("1. Specifying a model")
# Two factors with persistence 0.98 (slow) and 0.70 (fast), carrying 30% and
# 40% of the variance; the remaining 30% is white noise.
g <- arfactcop(lambda = c(0.98, 0.70), w = c(0.30, 0.40))
g
coef(g)
# The autocorrelations are rho_k = 0.3 * 0.98^k + 0.4 * 0.70^k (k >= 1):
round(0.3 * 0.98^(1:6) + 0.4 * 0.70^(1:6), 3)
# The same Gaussian copula as an ARMA(2,2) model: the AR polynomial is
# (1 - 0.98 B)(1 - 0.70 B) and the MA part comes from a spectral factorisation.
a <- arfact2arma(g)
a
# identical partial autocorrelations (and therefore kpacf and likelihood)
max(abs(ktau(g, 30) - ktau(a, 30)))
# A single factor plus noise is the ARMA(1,1) of GARCH(1,1) squares: with
# alpha1 = 0.1 and beta1 = 0.85 the persistence is 0.95 and the signal share 0.1885,
# and the implied ARMA(1,1) has ar = 0.95 and ma = -beta1 = -0.85:
arfact2arma(arfactcop(lambda = 0.95, w = 0.1884945))
# Structures to be fitted, and a model without a noise component (weights sum to one):
arfactcop(order = 2)
arfactcop(lambda = c(0.9, 0.5), w = c(0.6, 0.4), noise = FALSE)
try(arfactcop(lambda = c(0.9, 1.2), w = c(0.3, 0.3)))

# ---------------------------------------------------------------------------
section("2. The kpacf: fast decay first, slow decay later")
# Kendall's tau of the pair copulas of the equivalent D-vine, (2/pi) asin of the
# partial autocorrelations: non-negative at every lag.
tau <- ktau(g, 40)
round(tau[1:8], 3)
one <- arfactcop(lambda = 0.98, w = 0.7)               # one factor
barplot(rbind(ktau(one, 25), tau[1:25]), beside = TRUE, names.arg = 1:25,
        xlab = "lag", ylab = "Kendall's tau", legend.text = c("one factor", "two factors"),
        main = "kpacf of the factor models")

# ---------------------------------------------------------------------------
section("3. Simulation (tssim) and the exact Gaussian likelihood")
u <- tssim(g, 1500)
par(mfrow = c(1, 3))
plot(u, type = "l", main = "simulated copula data", xlab = "t", ylab = "u")
acf(qnorm(u), main = "ACF of the normal scores")
plot(0:30, c(1, 0.3 * 0.98^(1:30) + 0.4 * 0.7^(1:30)), type = "h", main = "model ACF", xlab = "lag", ylab = "")
par(op)

# ---------------------------------------------------------------------------
section("4. Fitting (tsfit)")
# Starting values come from a least-squares fit to the sample autocorrelations
# of the normal scores; the likelihood is exact (Kalman filter) and handles NA.
fit2 <- tsfit(u, arfactcop(order = 2), pseudo = FALSE)
fit2
arfact2arma(fit2)     # the implied ARMA(2,2) coefficients of the fit
fit1 <- tsfit(u, arfactcop(order = 1), pseudo = FALSE)
c(one_factor = AIC(fit1), two_factors = AIC(fit2))
# the best ARMA(1,1) for comparison (one factor plus noise is a special ARMA(1,1))
c(armacop_1_1 = AIC(tsfit(u, armacop(order = c(1, 1)), pseudo = FALSE)))

# ---------------------------------------------------------------------------
section("5. Standard errors, residuals and prediction")
fit_h <- tsfit(u, arfactcop(order = 2), pseudo = FALSE, se = "hessian")
fit_h
# Residuals (standardised one-step prediction errors) are iid under the model
r <- resid(fit2, scale = "normal")
par(mfrow = c(1, 2)); acf(r, main = "ACF of residuals"); qqnorm(r); qqline(r); par(op)
# Conditional mean and sd of the normal score h steps ahead, and the one-step predictive quantiles
predict(fit2, h = 1:5)
predict(fit2, type = "qf", x = c(0.05, 0.5, 0.95))

# ---------------------------------------------------------------------------
section("6. As the kpacf of a d-vine with absolute spherical t copulas")
# The partial autocorrelations of the factor model can calibrate the pair
# copulas of a d-vine (here ast copulas at lags 1 and 2, Gaussian beyond).
m <- sdvinecop(g, pairs = list(pc("ast"), pc("ast")))
m
ud <- tssim(m, 600)
fd <- tsfit(ud, sdvinecop(arfactcop(order = 2), pairs = list(pc("ast"), pc("ast"))), pseudo = FALSE)
fd
# automatic selection of the pair copulas, starting from a fitted factor model
auto <- autodvinecop(ud, tsfit(ud, arfactcop(order = 2), pseudo = FALSE), nreplace = 3,
                     choices = c("gumbel", "joe", "ast"), pseudo = FALSE)
cat("\nEnd of demo.\n")
