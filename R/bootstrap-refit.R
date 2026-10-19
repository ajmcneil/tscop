# Fast bootstrap refits for the d-vine and v-transform copula models.
#
# A bootstrap series is simulated from the fitted model and refitted from the
# fitted values, so the replicate's optimum lies close to the starting point and
# the curvature there is that of the fitted model. The refit is a BFGS iteration
# whose initial inverse Hessian is the inverse Hessian of the objective at the
# original fit (computed once from central differences of the exact gradient)
# instead of the identity: fewer iterations and, above all, far fewer line-search
# value evaluations than optim()'s BFGS, which has to rebuild the curvature.
# Replicates for which the iteration does not converge (or if the Hessian is not
# positive definite, for example at a boundary or a cusp of the likelihood) are
# refitted by BFGS as before. options(tscop.bootnewton = FALSE) forces BFGS.

# Hessian of the objective at z by central differences of the gradient `gr`
# (symmetrised); may contain non-finite values.
num_hessian_gr <- function(gr, z, h = 1e-4) {
  nz <- length(z)
  H <- matrix(0, nz, nz)
  for (k in seq_len(nz)) {
    zp <- z; zp[k] <- zp[k] + h
    zm <- z; zm[k] <- zm[k] - h
    H[, k] <- (gr(zp) - gr(zm)) / (2 * h)
  }
  (H + t(H)) / 2
}

# That Hessian, or NULL when it is not positive definite or the gradient is not finite.
boot_hessian <- function(gr, z, h = 1e-4) {
  H <- num_hessian_gr(gr, z, h)
  if (any(!is.finite(H))) return(NULL)
  ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) <= 1e-8 * max(ev)) return(NULL)
  H
}

# Inverse Hessian for the Newton refits, or NULL (use BFGS).
boot_prepare <- function(gr, z) {
  if (!isTRUE(getOption("tscop.bootnewton", TRUE)) || length(z) == 0L) return(NULL)
  H <- tryCatch(boot_hessian(gr, z), error = function(e) NULL)
  if (is.null(H)) NULL else solve(H)
}

# Quasi-Newton (BFGS) minimisation of fn from z0 whose initial inverse Hessian is Hinv,
# the inverse Hessian at the original fit, with the exact gradient gr and a backtracking
# line search. With the curvature of the fitted model as the starting point the first
# step is nearly always accepted, and the updates adapt the matrix to the replicate.
# Stops when the Newton decrement (the predicted decrease of fn) falls below tol; returns
# list(par, value), or NULL if it does not converge.
boot_newton <- function(z0, fn, gr, Hinv, tol = 1e-8, maxit = 40L) {
  z <- z0
  f <- fn(z)
  if (!is.finite(f) || f >= 1e10) return(NULL)
  g <- gr(z)
  if (any(!is.finite(g))) return(NULL)
  Hk <- Hinv
  nz <- length(z)
  for (it in seq_len(maxit)) {
    step <- -drop(Hk %*% g)
    dec <- -sum(g * step)
    if (!is.finite(dec) || dec <= 0) {          # not a descent direction: restart from the fitted curvature
      Hk <- Hinv
      step <- -drop(Hk %*% g)
      dec <- -sum(g * step)
    }
    if (dec < 2 * tol) return(list(par = z, value = f))
    tt <- 1
    repeat {
      zn <- z + tt * step
      fn_new <- fn(zn)
      if (is.finite(fn_new) && fn_new < f - 1e-4 * tt * dec) break
      tt <- tt / 2
      if (tt < 1e-4) return(NULL)
    }
    gn <- gr(zn)
    if (any(!is.finite(gn))) return(NULL)
    s <- zn - z
    y <- gn - g
    sy <- sum(s * y)
    if (sy > 1e-10 * sqrt(sum(s * s) * sum(y * y))) {
      rho <- 1 / sy
      I <- diag(nz)
      Hk <- (I - rho * tcrossprod(s, y)) %*% Hk %*% (I - rho * tcrossprod(y, s)) + rho * tcrossprod(s)
    }
    z <- zn
    f <- fn_new
    g <- gn
  }
  NULL
}

# Run the B replicates, in parallel when cores > 1 (not on Windows); `refit(b)` returns
# the coefficient vector or NULL. The replicates use no random numbers, so the result does
# not depend on `cores`.
boot_run <- function(B, refit, cores = 1L) {
  if (cores > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(seq_len(B), refit, mc.cores = cores)
  } else {
    lapply(seq_len(B), refit)
  }
}
