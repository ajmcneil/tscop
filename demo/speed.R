# Compiled kernels: what they cover and what they buy. The likelihood of the d-vine
# models (sdvinecop, vtscop with a d-vine base), its exact gradient, the simulation, the
# Rosenblatt transform and the one-step prediction run in compiled code for the pair
# copulas Gauss, Clayton, Gumbel, Joe, Frank, BB1, t and the absolute spherical t (ast), with
# rotations, and for the inverse-v-transformed (IVT) pairs with the bases ast, joe, gumbel,
# clayton180, frank and gauss; the Kalman filters of armacop and arfactcop are compiled too.
# Everything else (a rotated ast pair) falls back to the R implementation, which is kept as
# the reference and can be selected with options(tscop.cpp = FALSE).
# Run with demo("speed", package = "tscop") or
#   devtools::load_all(); source("demo/speed.R", echo = TRUE)
# Takes about 25 seconds (most of it the R implementation, for comparison).

if (!("package:tscop" %in% search())) library(tscop)
section <- function(title) cat("\n\n==== ", title, " ====\n\n", sep = "")
elapsed <- function(expr) unname(system.time(expr)[3])
ncores <- if (.Platform$OS.type == "windows") 1L else max(1L, parallel::detectCores() - 1L)
set.seed(11)

# ---------------------------------------------------------------------------
section("1. A fit with and without the compiled kernel")
truth <- sdvinecop(arfactcop(lambda = 0.95, w = 0.5), tail = pc("ivt", base = "joe"),
                   vt = c(0.45, 0.55), maxlag = 40, tautol = 1e-3)
u <- tssim(truth, 2000)
spec <- sdvinecop(arfactcop(order = 1), tail = pc("ivt", base = "joe"), maxlag = 40, tautol = 1e-3)
t_cpp <- elapsed(f_cpp <- tsfit(u, spec, pseudo = FALSE))
options(tscop.cpp = FALSE)                       # the R implementation
t_r <- elapsed(f_r <- tsfit(u, spec, pseudo = FALSE))
options(tscop.cpp = TRUE)
c(compiled = t_cpp, R = t_r, speedup = t_r / t_cpp)
# the same optimum
rbind(compiled = coef(f_cpp), R = coef(f_r))
c(compiled = as.numeric(logLik(f_cpp)), R = as.numeric(logLik(f_r)))

# ---------------------------------------------------------------------------
section("2. The exact gradient")
# BFGS uses the exact gradient from forward-mode derivatives in the kernel; options(tscop.cppgrad
# = FALSE) uses central differences of the compiled likelihood instead.
options(tscop.cppgrad = FALSE)
t_fd <- elapsed(f_fd <- tsfit(u, spec, pseudo = FALSE))
options(tscop.cppgrad = TRUE)
c(exact_gradient = t_cpp, central_differences = t_fd)

# ---------------------------------------------------------------------------
section("3. Simulation and the absolute spherical t copula")
# The ast pair copula needs a Student t quantile function at every pair and level; the kernel
# carries t-scale values through tabulated maps (see ?pc and the notes in docs/).
ast <- sdvinecop(arfactcop(lambda = 0.95, w = 0.5), tail = pc("ivt", base = "ast"),
                 vt = c(0.45, 0.55), maxlag = 40, tautol = 1e-3)
t_sim <- elapsed(x_cpp <- tssim(ast, 3000))
options(tscop.cpp = FALSE)
t_sim_r <- elapsed(x_r <- tssim(ast, 1000))      # (a shorter series: the R version is slow)
options(tscop.cpp = TRUE)
cat(sprintf("simulating an IVT-ast series: %.2f s for 3000 values (compiled); %.1f s for 1000 values (R)\n", t_sim, t_sim_r))
# fitting it takes about a second compiled (the R implementation needs about 25 s for this fit)
t_ast <- elapsed(f_ast <- tsfit(x_cpp, sdvinecop(arfactcop(order = 1), tail = pc("ivt", base = "ast"), maxlag = 40, tautol = 1e-3),
                               pseudo = FALSE))
cat(sprintf("fitting the IVT-ast model to 3000 values: %.1f s\n", t_ast))
coef(f_ast)

# ---------------------------------------------------------------------------
section("4. Bootstrap standard errors on several cores")
# The B refits start from the fitted values with the curvature of the fit as the initial
# inverse Hessian; they use no random numbers, so the results do not depend on cores. (The gain
# from several cores is below the number of cores: the replicates are simulated and the
# Hessian computed in the main process, and forked processes have a start-up cost.)
B <- 100
t1 <- elapsed(b1 <- tsfit(u, spec, pseudo = FALSE, se = "bootstrap", B = B))
tn <- elapsed(bn <- tsfit(u, spec, pseudo = FALSE, se = "bootstrap", B = B, cores = ncores))
c(one_core = t1, several_cores = tn, cores = ncores)
rbind(se_one_core = b1@se, se_several_cores = bn@se)

# ---------------------------------------------------------------------------
section("5. Full time series models: joint maximum likelihood")
tm <- tsmod(marg("st", c(mu = 0, sigma = 1, df = 5)),
            sdvinecop(arfactcop(lambda = 0.95, w = 0.5), pairs = list(pc("bb1", aux = 1.2)), tail = pc("joe"),
                      maxlag = 30, tautol = 1e-3))
x <- tssim(tm, 3000)
tspec <- tsmod("st", sdvinecop(arfactcop(order = 1), pairs = list(pc("bb1", aux = 1.1)), tail = pc("joe"),
                              maxlag = 30, tautol = 1e-3))
t_two <- elapsed(f_two <- tsfit(x, tspec))
t_joint <- elapsed(f_joint <- tsfit(x, tspec, method = "joint"))
c(two_stage = t_two, joint_including_the_two_stage_start = t_joint)
c(two_stage = as.numeric(logLik(f_two)), joint = as.numeric(logLik(f_joint)))
rbind(two_stage = coef(f_two), joint = coef(f_joint))
cat("\nEnd of demo.\n")
