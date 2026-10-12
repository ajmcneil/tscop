# Parametric margins (marg): Gaussian, Laplace, skew Laplace, Student t and
# skew Student t, with sigma a pure scale parameter and the skew families built
# by Fernandez-Steel skewing (gamma = 1 symmetric, mu the mode). Simulation,
# fitting by maximum likelihood with standard errors, fixed parameters,
# residuals (probability integral transforms) and comparison of nested families.
# Run with demo("marg", package = "tscop") or
#   devtools::load_all(); source("demo/marg.R", echo = TRUE)

if (!("package:tscop" %in% search())) library(tscop)
section <- function(title) cat("\n\n==== ", title, " ====\n\n", sep = "")
op <- par(no.readonly = TRUE)
set.seed(1)

section("1. Specification, functions and simulation")
m <- marg("sst", c(mu = 0, sigma = 1, df = 5, gamma = 1.3))
m
mdens(m, c(-2, 0, 2))                       # density
mdf(m, c(-2, 0, 2))                         # distribution function
mqf(m, c(0.01, 0.5, 0.99))                  # quantile function
x <- tssim(m, 4000)
par(mfrow = c(1, 2))
hist(x, breaks = 60, freq = FALSE, main = "simulated sst data", xlab = "x")
curve(mdens(m, x), add = TRUE, col = 2, lwd = 2)
plot(density(x), main = "log density", log = "y", xlab = "x")
curve(mdens(m, x), add = TRUE, col = 2)
par(op)

section("2. Fitting with standard errors")
fit <- tsfit(x, marg("sst"), se = "hessian")
fit
rbind(true = m@pars, estimate = fit@model@pars)
fitb <- tsfit(x, marg("sst"), se = "bootstrap", B = 50)   # parametric bootstrap
round(cbind(hessian = fit@se, bootstrap = fitb@se), 4)

section("3. Fixed parameters")
# a zero-mean Student t: mu is held at zero and not estimated
zt <- tsfit(x, marg("st", c(mu = 0), fix = "mu"))
zt
c(npar = attr(logLik(zt), "df"))

section("4. Comparing nested families")
fits <- list(gauss = tsfit(x, marg("gauss")), laplace = tsfit(x, marg("laplace")),
             slaplace = tsfit(x, marg("slaplace")), st = tsfit(x, marg("st")),
             sst = fit)
aic <- sapply(fits, AIC)
round(aic - min(aic), 1)                      # sst, the true family, should be best

section("5. Residuals")
r <- resid(fit)                               # probability integral transforms
par(mfrow = c(1, 2))
hist(r, main = "PITs of the fitted sst", xlab = "u")
qqnorm(resid(fit, scale = "normal")); abline(0, 1)
par(op)
ks.test(r, "punif")$p.value
cat("\nEnd of demo.\n")
