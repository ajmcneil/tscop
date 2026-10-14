# Demonstration of the sdvinecop class: stationary d-vine copula processes with
# an ARMA kpacf and non-Gaussian pair copulas. Run it with
#   demo("sdvinecop", package = "tscop")      (after installing tscop), or
#   devtools::load_all(); source("demo/sdvinecop.R", echo = TRUE)
# from the source folder. It produces several plots and takes about twenty seconds
# (the likelihood, its gradient, the simulation and the bootstrap run in compiled code;
# see demo("speed", package = "tscop")).

if (!("package:tscop" %in% search())) library(tscop)
set.seed(2027)
ncores <- if (.Platform$OS.type == "windows") 1L else max(1L, parallel::detectCores() - 1L)   # processes for the bootstrap refits
op <- par(no.readonly = TRUE)

section <- function(title) cat("\n\n==== ", title, " ====\n\n", sep = "")

# ---------------------------------------------------------------------------
section("1. Specifying a model")
# A tied d-vine: a Gaussian ARMA model supplies the Kendall partial
# autocorrelations tau_k (the kpacf); the pair copula at lag k is a member of a
# chosen family calibrated to have Kendall's tau tau_k. Lags without an
# explicit specification use the tail family (Gaussian by default).
arma <- armacop(ar = 0.5, ma = 0.4)
mod <- sdvinecop(arma, pairs = list(
  pc("gumbel"),                           # lag 1: upper-tail dependent
  pc("clayton", pos = 180, neg = 270),    # lag 2: survival Clayton / Clayton 270
  pc("t", aux = 5)))                      # lag 3: t copula, 5 degrees of freedom
mod
# pc() stores the rotations for positive (pos) and negative (neg) dependence;
# the one used at each lag follows the sign of its tau. Both cases are shown by
# show(): lag 2 here has a negative tau, so Clayton 270 is used.
coef(mod)

# The effective lag K is the last lag at which |tau_k| exceeds tautol; the
# model is a finite-order Markov chain of order K. It is stored in the object.
mod@K
sdvinecop(arma, pairs = list(pc("gumbel")), tautol = 1e-2)@K     # coarser truncation
sdvinecop(arma, pairs = list(pc("gumbel")), maxlag = 5)@K        # imposed cap
# A pure AR kpacf has exactly the AR order:
sdvinecop(armacop(ar = c(0.5, 0.2)), pairs = list(pc("gumbel"), pc("frank")))@K

# Tail families that tend to independence are allowed as the tail:
sdvinecop(arma, pairs = list(pc("gumbel")), tail = pc("clayton", pos = 180))
try(sdvinecop(arma, tail = "t"))   # t and BB1 cannot be tails

# ---------------------------------------------------------------------------
section("2. Kendall's tau of the pair copulas (ktau)")
tau <- ktau(mod, 20)
round(tau, 3)
barplot(tau, names.arg = 1:20, xlab = "lag", ylab = "Kendall's tau",
        main = "kpacf of the d-vine (zero beyond K)")

# ---------------------------------------------------------------------------
section("3. Simulation (tssim) and non-Gaussian features")
# Exact stationary simulation (inverse Rosenblatt transform, O(nK)).
u <- tssim(mod, 600)
g <- tssim(sdvinecop(arma), 600)    # the Gaussian ARMA copula with the same kpacf
par(mfrow = c(2, 2))
plot(u, type = "l", xlab = "t", ylab = "u", main = "d-vine with Gumbel/Clayton/t")
plot(g, type = "l", xlab = "t", ylab = "u", main = "Gaussian ARMA copula")
plot(u[-600], u[-1], pch = 20, cex = 0.5, xlab = "u[t]", ylab = "u[t+1]",
     main = "lag-1 scatter: upper-tail dependence")
plot(g[-600], g[-1], pch = 20, cex = 0.5, xlab = "u[t]", ylab = "u[t+1]",
     main = "lag-1 scatter: Gaussian")
par(op)

# ---------------------------------------------------------------------------
section("4. Residuals (resid) and prediction (predict)")
# Rosenblatt residuals are iid uniform from the first observation if the model
# is right; compare the true and a Gaussian model.
r_true <- resid(mod, u, scale = "normal")
r_gauss <- resid(sdvinecop(arma), u, scale = "normal")
par(mfrow = c(1, 3))
qqnorm(r_true, main = "true model"); qqline(r_true)
qqnorm(r_gauss, main = "Gaussian model"); qqline(r_gauss)
acf(r_true, main = "ACF of residuals (true model)")
par(op)
# One-step predictive distribution after a large value: tail dependence makes
# the Gumbel model put more mass in the upper tail than the Gaussian model.
x <- seq(0.01, 0.99, length.out = 199)
hist_hi <- c(u[1:(length(u) - 1)], 0.97)
plot(x, predict(mod, hist_hi, type = "dens", x = x), type = "l", lwd = 2,
     xlab = "u[t+1]", ylab = "predictive density",
     main = "one-step density after u[t] = 0.97")
lines(x, predict(sdvinecop(arma), hist_hi, type = "dens", x = x), lty = 2, col = 2)
legend("topleft", c("d-vine", "Gaussian"), lty = c(1, 2), col = c(1, 2), bty = "n")
predict(mod, u, type = "qf", x = c(0.05, 0.5, 0.95))   # predictive quantiles

# ---------------------------------------------------------------------------
section("5. Fitting a specified model (tsfit)")
u <- tssim(mod, 3000)                         # a series of realistic length
spec <- sdvinecop(armacop(order = c(1, 1)),
                 pairs = list(pc("gumbel"), pc("clayton", pos = 180, neg = 270), pc("t", aux = 5)))
t0 <- Sys.time()
fit <- tsfit(u, spec, pseudo = FALSE)
fit
cat("fitted in", format(round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)), "seconds\n")
# The effective lag was fixed at its starting value during the optimisation and
# updated by a rule at the end (refits recorded in the object):
c(Kstart = fit@Kstart, K = fit@model@K, refits = fit@nrefit)
# Compare with the Gaussian ARMA copula (AIC counts the t degrees of freedom)
gfit <- tsfit(u, sdvinecop(armacop(order = c(1, 1))), pseudo = FALSE)
c(dvine = AIC(fit), gaussian = AIC(gfit))
AICc(fit, gfit)                # AIC, BIC and the small-sample AICc are available for every fitted class

# ---------------------------------------------------------------------------
section("6. Standard errors")
fit_h <- tsfit(u, spec, pseudo = FALSE, se = "hessian")
fit_h
# Parametric bootstrap with B = 200 replicates, simulated together and refitted from
# the fitted values with K fixed (the refits run on several cores; the results do not
# depend on the number of cores).
system.time(fit_b <- tsfit(u, spec, pseudo = FALSE, se = "bootstrap", B = 200, cores = ncores))
rbind(hessian = fit_h@se, bootstrap = fit_b@se)
# The bootstrap is recommended: Hessian standard errors are unreliable where a
# rotated copula switches rotation (tau near 0), at the boundary of a BB1
# copula's feasible region, and for pseudo-observations.

# ---------------------------------------------------------------------------
section("6b. Diagnostics: empirical Kendall partial autocorrelations (ktauhat)")
# The data's partial dependence at each lag, conditioning on the in-between values through the
# model's own pair copulas, against the model's kpacf. The band is the 95% interval of the same
# statistic for B series simulated from the model (parameter uncertainty is not included).
kh <- ktauhat(fit, lagmax = 12, B = 200, cores = ncores)
round(kh, 3)
matplot(kh$lag, kh[, c("lo", "hi")], type = "l", lty = 2, col = "gray50", ylim = range(kh[, -c(1, 4)]),
        xlab = "lag", ylab = "Kendall's tau", main = "empirical vs model kpacf")
points(kh$lag, kh$tau_data, pch = 19)
lines(kh$lag, kh$tau_model, col = 2, lwd = 2)
# A misspecified model is flagged: a smaller AR coefficient than the data's
bad <- sdvinecop(armacop(ar = 0.2, ma = 0.4), pairs = list(pc("gumbel"), pc("clayton", pos = 180, neg = 270), pc("t", aux = 5)))
round(ktauhat(bad, u, lagmax = 6, B = 100, cores = ncores), 3)

# ---------------------------------------------------------------------------
section("7. Automatic selection of the pair copulas (autodvinecop)")
# Start from the Gaussian ARMA copula and let the search choose copulas lag by
# lag (greedy, by AIC). Only the rotations for the realised sign of tau are
# tried; candidates are screened by one likelihood evaluation each, and the
# best two are refitted jointly. (Gumbel and survival Clayton copulas, like
# Gauss and Frank, are easily confused, especially at small |tau|; the aim is
# to estimate the kpacf well rather than to identify every family.)
gauss <- tsfit(u, armacop(order = c(1, 1)), pseudo = FALSE)
t0 <- Sys.time()
auto <- autodvinecop(u, gauss, nreplace = 8, aux = "t", pseudo = FALSE)
cat("search took", format(round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)), "seconds\n")
auto@path
c(gaussian = AIC(gfit), selected = AIC(auto), true_structure = AIC(fit))
# Screening off reproduces the original full-refit behaviour (slower):
system.time(autodvinecop(u, gauss, nreplace = 3, screen = FALSE, pseudo = FALSE, verbose = FALSE))
system.time(autodvinecop(u, gauss, nreplace = 3, screen = TRUE, pseudo = FALSE, verbose = FALSE))
# A search can be continued:
more <- autodvinecop(u, auto, nreplace = 7, pseudo = FALSE, verbose = FALSE)
more@path

# ---------------------------------------------------------------------------
section("8. Free mode: explicit pair copulas for simulation")
# A free-mode model lists rvinecopulib pair copulas with their own parameters
# (optionally followed by a Gaussian ARMA tail). It can be simulated from and
# evaluated, but not fitted: whether a freely specified set of one-signed
# copulas can match the pattern of positive and negative Kendall partial
# autocorrelations of data is not known in advance, so fitting uses tied mode.
free <- sdvinecop(pairs = list(rvinecopulib::bicop_dist("gumbel", 0, 2),
                              rvinecopulib::bicop_dist("clayton", 270, 1.5),
                              rvinecopulib::bicop_dist("frank", 0, 3)))
free
round(ktau(free, 5), 3)
uf <- tssim(free, 400)
plot(uf[-400], uf[-1], pch = 20, cex = 0.5, main = "free-mode d-vine, lag-1 scatter")
try(tsfit(uf, free))

# ---------------------------------------------------------------------------
section("9. Time series input and warm starts")
uts <- ts(u, start = c(2010, 1), frequency = 12)
fit_ts <- tsfit(uts, spec, pseudo = FALSE)
head(resid(fit_ts))          # the time index is carried to the residuals
# Refit the same structure to new data, starting from the previous estimates:
system.time(tsfit(tssim(fit@model, 3000), fit@model, pseudo = FALSE))
# Simulate from the fitted model
head(tssim(fit, 5))

# ---------------------------------------------------------------------------
section("10. Pseudo-observations")
x <- qgamma(u, shape = 2)                 # a non-normal margin
pobs <- rank(x) / (length(x) + 1)
tsfit(pobs, spec, pseudo = TRUE)
cat("\nEnd of demo.\n")
