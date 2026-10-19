# Joint maximum likelihood for tsmod models: the margin and the copula process
# are estimated together by maximising
#   sum_t log f(x_t; theta_m) + log c(F(x_t; theta_m); theta_c)
# on the free scale of both parameter sets. The copula enters only through the
# internal cop_obj() interface (R/copula-objective.R). The effective lag K of a
# d-vine copula is held fixed at its two-stage value; the number of nonlinear lags L is held fixed during an
# optimisation and follows the rule of the two-stage fit afterwards (tm_joint_L_rule). Margin parameters that put a PIT within 1e-9 of 0 or 1,
# or the median PIT outside (0.05, 0.95), are infeasible (the copula likelihood is unbounded when the PITs
# collapse onto the boundary).

# Joint objective for the fitted margin `mmodel` and copula `cmodel`.
tm_joint_setup <- function(mmodel, cmodel, x) {
  fam <- mg_family(mmodel@family)
  th0 <- mmodel@pars
  free <- setdiff(fam$pars, mmodel@fix)
  types <- fam$types[free]
  nm <- length(free)
  u0 <- tm_pits(mmodel, x)
  cop <- cop_obj(cmodel, u0)
  if (!is.null(cop$setL) && is.na(cop$L())) cop$setL(cop$Lat(cop$z0))     # L is held fixed during the optimisation
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
    if (any(u < 1e-9 | u > 1 - 1e-9) || stats::median(u) < 0.05 || stats::median(u) > 0.95) return(1e10)
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

# Best of the points par + s e_k (coordinate steps s = +-1e-4, ..., +-0.1 on the free scale) and `par`.
tm_joint_probe <- function(jt, x, par, value) {
  best <- list(par = par, value = value)
  for (k in seq_along(par)) {
    for (s in c(-1, 1) * rep(10^(-4:-1), each = 2)) {
      p <- par; p[k] <- p[k] + s
      v <- jt$negll(p, x)
      if (is.finite(v) && v < best$value) best <- list(par = p, value = v)
    }
  }
  best
}

# Maximise the joint likelihood from `start`. The result never has a lower likelihood than the start. The objective can
# be rough in the margin parameters (observations cross the fulcrums of an inverse v-transform as the PITs shift), so that
# BFGS can stall after a tiny step and still report convergence. After BFGS the point is therefore probed with coordinate
# steps of several sizes; if one of them improves the log-likelihood by `tol` (0.01) or more, BFGS is restarted from there
# on a coarser parameter scale (at most three times), and Nelder-Mead is the last resort. `convergence` is 0 (converged:
# no probe improves by `tol`), 1 (iteration limit) or 10 (stalled: the probes still find improvements).
tm_joint_optim <- function(jt, x, start = jt$theta0, control = list(), tol = 0.01) {
  np <- length(start)
  v0 <- jt$negll(start, x)
  if (np == 0L) return(list(theta = start, value = v0, convergence = 0L, message = ""))
  if (v0 >= 1e10) stop("the two-stage estimates are infeasible for the joint likelihood", call. = FALSE)
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
  best <- list(par = start, value = v0)
  conv <- 0L
  msg <- ""
  run <- function(from, ct, method = "BFGS") {
    opt <- if (method == "BFGS") {
      stats::optim(from, jt$negll, gr = gr, x = x, method = "BFGS", control = ct)
    } else {
      stats::optim(from, jt$negll, x = x, method = "Nelder-Mead", control = ct)
    }
    if (is.finite(opt$value) && opt$value < best$value) best <<- list(par = opt$par, value = opt$value)
    conv <<- opt$convergence
    msg <<- if (is.null(opt$message)) "" else opt$message
  }
  c2 <- utils::modifyList(ctrl, list(parscale = rep(0.05, np)))
  run(start, ctrl)
  ok <- FALSE
  for (it in 1:4) {
    pr <- tm_joint_probe(jt, x, best$par, best$value)
    if (best$value - pr$value < tol) {
      ok <- TRUE
      break
    }
    if (it == 4L) break
    best <- pr
    run(best$par, c2)
  }
  if (!ok || conv != 0L) {
    run(best$par, list(maxit = 2000, reltol = 1e-10, warn.1d.NelderMead = FALSE), "NM")
    pr <- tm_joint_probe(jt, x, best$par, best$value)
    if (best$value - pr$value >= tol) {
      best <- pr
      conv <- 10L
      msg <- "the optimiser stalled: a coordinate step still improves the log-likelihood"
    }
  }
  list(theta = best$par, value = best$value, convergence = as.integer(conv), message = msg)
}

# The number of nonlinear lags L (nonlinear v-transforms) at the joint optimum follows the rule of the two-stage fit: if
# the rule gives another number at the optimum, the log-likelihood is compared with the new one at the same parameters; a
# change below `tol` keeps the new number, otherwise the fit is repeated with it from the optimum and kept only if it has
# the higher likelihood (at most twice).
tm_joint_L_rule <- function(jt, x, fit, control = list(), tol = 0.01) {
  cop <- jt$cop
  if (is.null(cop$Lat) || is.null(cop$setL) || is.na(cop$L())) return(fit)
  zc <- function(th) th[jt$nm + seq_along(cop$z0)]
  for (it in 1:2) {
    Lold <- cop$L()
    Lnew <- cop$Lat(zc(fit$theta))
    if (is.na(Lnew) || Lnew == Lold) break
    ll_old <- -jt$negll(fit$theta, x)
    cop$setL(Lnew)
    ll_new <- -jt$negll(fit$theta, x)
    if (is.finite(ll_new) && abs(ll_new - ll_old) < tol) {
      fit$value <- -ll_new
      break
    }
    new <- tryCatch(tm_joint_optim(jt, x, start = fit$theta, control = control), error = function(e) NULL)
    if (is.null(new) || new$value >= fit$value) {
      cop$setL(Lold)
      break
    }
    fit <- new
  }
  fit
}

# Full joint fit: the start is the two-stage fit (`cfit_model`, the copula refitted to the PITs of `mmodel`), or the
# copula model supplied by the user (`cmodel`) with the same margin if that has the higher likelihood (a fitted model is
# then a start as good as the one it came from, instead of a refit that can end in a worse local optimum);
# the optimisation, then the rule for L. Returns the setup used and the fit, with `start` naming the start.
tm_joint_fit <- function(mmodel, cmodel, cfit_model, x, control = list()) {
  jt <- tm_joint_setup(mmodel, cfit_model, x)
  start <- "twostage"
  if (!identical(cmodel, cfit_model)) {
    alt <- tryCatch(tm_joint_setup(mmodel, cmodel, x), error = function(e) NULL)
    if (!is.null(alt) && length(alt$theta0) == length(jt$theta0)) {
      va <- alt$negll(alt$theta0, x)
      if (is.finite(va) && va < 1e10 && va < jt$negll(jt$theta0, x) - 1e-8) {
        jt <- alt
        start <- "supplied"
      }
    }
  }
  fit <- tm_joint_optim(jt, x, control = control)
  fit <- tm_joint_L_rule(jt, x, fit, control)
  c(fit, list(jt = jt, start = start))
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
