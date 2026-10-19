# Joint maximum likelihood for tsmod models: the margin and the copula process
# are estimated together by maximising
#   sum_t log f(x_t; theta_m) + log c(F(x_t; theta_m); theta_c)
# on the free scale of both parameter sets. The copula enters only through the
# internal cop_obj() interface (R/copula-objective.R). The effective lag of a
# d-vine copula is held fixed at its two-stage value. Margin parameters that put a PIT within 1e-9 of 0 or 1
# are infeasible (the copula likelihood is unbounded when the PITs collapse onto the boundary).

# Joint objective for the fitted margin `mmodel` and copula `cmodel`.
tm_joint_setup <- function(mmodel, cmodel, x) {
  fam <- mg_family(mmodel@family)
  th0 <- mmodel@pars
  free <- setdiff(fam$pars, mmodel@fix)
  types <- fam$types[free]
  nm <- length(free)
  u0 <- tm_pits(mmodel, x)
  cop <- cop_obj(cmodel, u0)
  zm0 <- if (nm > 0L) mg_to_free(th0[free], types) else numeric(0)
  mpars <- function(zm) {
    th <- th0
    if (nm > 0L) th[free] <- mg_from_free(zm, types)
    th
  }
  negll <- function(theta, x) {
    zm <- theta[seq_len(nm)]
    zc <- theta[nm + seq_along(cop$z0)]
    vm <- if (nm > 0L) mg_negll(zm, fam, free, th0, x) else -sum(mg_ldens(fam, th0, x))
    if (!is.finite(vm) || vm >= 1e10) return(1e10)
    u <- mg_cdf(fam, mpars(zm), x)
    if (anyNA(u)) return(1e10)
    # a margin that pushes observations to the ends of (0, 1) makes the copula likelihood unbounded
    # (all the PITs coincide at the clamp): such parameter values are infeasible
    if (any(u < 1e-9 | u > 1 - 1e-9)) return(1e10)
    u <- pmin(pmax(u, 1e-12), 1 - 1e-12)
    vc <- tryCatch(cop$negll(zc, u), error = function(e) 1e10)
    if (!is.finite(vc) || vc >= 1e10) return(1e10)
    vm + vc
  }
  # Exact gradient of the joint negative log-likelihood where the copula supplies one: the
  # margin part by central differences of its density (cheap), the copula part by the kernel with
  # the derivatives of the PITs u = F(x; theta_m) with respect to the free margin parameters
  # passed as derivatives of the data. NULL (differences are used) otherwise.
  gradient <- function(theta, x) {
    if (is.null(cop$grad) || !isTRUE(getOption("tscop.cppgrad", TRUE))) return(NULL)
    zm <- theta[seq_len(nm)]
    zc <- theta[nm + seq_along(cop$z0)]
    if (any(!is.finite(theta))) return(NULL)
    clampu <- function(v) pmin(pmax(v, 1e-12), 1 - 1e-12)
    u <- clampu(mg_cdf(fam, mpars(zm), x))
    if (anyNA(u)) return(NULL)
    gm <- numeric(nm)
    du <- matrix(0, length(x), nm)
    h <- 1e-5
    for (k in seq_len(nm)) {
      zp <- zm; zp[k] <- zp[k] + h
      zq <- zm; zq[k] <- zq[k] - h
      gm[k] <- (mg_negll(zp, fam, free, th0, x) - mg_negll(zq, fam, free, th0, x)) / (2 * h)
      du[, k] <- (clampu(mg_cdf(fam, mpars(zp), x)) - clampu(mg_cdf(fam, mpars(zq), x))) / (2 * h)
    }
    if (nm > 0L && (any(!is.finite(gm)) || any(!is.finite(du)))) return(NULL)
    gc <- cop$grad(zc, u, du)
    if (is.null(gc) || any(!is.finite(gc))) return(NULL)
    c(gm, rep(0, length(cop$z0))) + gc
  }
  coefmap <- function(theta) {
    zm <- theta[seq_len(nm)]
    zc <- theta[nm + seq_along(cop$z0)]
    c(stats::setNames(if (nm > 0L) mg_from_free(zm, types) else numeric(0), if (nm > 0L) paste0("margin.", free)),
      stats::setNames(cop$coef(zc), paste0("copula.", names(cop$coef(zc)))))
  }
  list(fam = fam, th0 = th0, free = free, nm = nm, cop = cop, mmodel = mmodel,
       theta0 = c(zm0, cop$z0), negll = negll, gradient = gradient, coefmap = coefmap,
       build = function(theta) {
         zm <- theta[seq_len(nm)]
         list(margin = new("marg", family = mmodel@family, pars = mpars(zm), fix = mmodel@fix),
              copula = cop$build(theta[nm + seq_along(cop$z0)]))
       })
}

# Maximise the joint likelihood from `start` (BFGS, Nelder-Mead fallback).
tm_joint_optim <- function(jt, x, start = jt$theta0, control = list()) {
  np <- length(start)
  if (np == 0L) return(list(theta = start, value = jt$negll(start, x), convergence = 0L, message = ""))
  if (jt$negll(start, x) >= 1e10) stop("the two-stage estimates are infeasible for the joint likelihood", call. = FALSE)
  ctrl <- utils::modifyList(list(maxit = 500, reltol = 1e-10, ndeps = rep(1e-6, np)), control)
  gr <- function(theta, x) {
    g <- jt$gradient(theta, x)
    if (!is.null(g)) return(g)
    vapply(seq_along(theta), function(k) {
      h <- 1e-6
      tp <- theta; tp[k] <- tp[k] + h
      tm <- theta; tm[k] <- tm[k] - h
      (jt$negll(tp, x) - jt$negll(tm, x)) / (2 * h)
    }, numeric(1))
  }
  opt <- stats::optim(start, jt$negll, gr = gr, x = x, method = "BFGS", control = ctrl)
  if (opt$convergence != 0L || opt$value >= 1e10) {
    alt <- stats::optim(opt$par, jt$negll, x = x, method = "Nelder-Mead",
                        control = list(maxit = 5000, reltol = 1e-10, warn.1d.NelderMead = FALSE))
    if (alt$value < opt$value) opt <- alt
  }
  list(theta = opt$par, value = opt$value, convergence = opt$convergence,
       message = if (is.null(opt$message)) "" else opt$message)
}

# A margfit object for the jointly estimated margin (standard errors are
# filled in by the caller through the combined covariance matrix).
tm_joint_margfit <- function(jt, margin, x, start_model, convergence, message) {
  free <- jt$free
  npar <- length(free)
  fam <- jt$fam
  new("margfit", model = margin, data = x, estimate = margin@pars[free],
      se = rep(NA_real_, npar), vcov = matrix(NA_real_, npar, npar, dimnames = list(free, free)),
      se_method = "none", boot = matrix(NA_real_, 0L, npar, dimnames = list(NULL, free)),
      loglik = sum(mg_ldens(fam, margin@pars, x)), nobs = length(x), npar = as.integer(npar),
      convergence = as.integer(convergence), message = message, start = start_model, template = NULL)
}
