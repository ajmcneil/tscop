# Reference implementation of the d-vine log-density for tests: a plain
# lag-by-lag loop using only rvinecopulib, with no normal-score shortcuts.
dv_loglik_ref <- function(u, pairs) {
  n <- length(u)
  K <- min(length(pairs), n - 1L)
  v <- cbind(u[-n], u[-1])
  ll <- 0
  for (j in seq_len(K)) {
    pr <- pairs[[j]]
    ll <- ll + sum(log(rvinecopulib::dbicop(v, pr$family, pr$rotation, pr$parameters)))
    if (j < K) {
      m <- nrow(v)
      v <- cbind(
        as.numeric(rvinecopulib::hbicop(v[-m, , drop = FALSE], 2L, pr$family, pr$rotation, pr$parameters)),
        as.numeric(rvinecopulib::hbicop(v[-1, , drop = FALSE], 1L, pr$family, pr$rotation, pr$parameters))
      )
    }
  }
  ll
}

# Log-density of u under the d-vine with the given resolved pair copulas, via
# rvinecopulib::vinecop_dist with a d-vine structure. rvinecopulib's
# structure 1:n pairs variables in time order and uses the same rotation
# convention as the engine.
dv_loglik_vinecop <- function(u, pairs) {
  n <- length(u)
  K <- length(pairs)
  pcs <- lapply(seq_len(n - 1L), function(j) {
    if (j <= K) {
      f <- pairs[[j]]
      replicate(n - j, rvinecopulib::bicop_dist(f$family, f$rotation, f$parameters),
                simplify = FALSE)
    } else {
      replicate(n - j, rvinecopulib::bicop_dist("indep"), simplify = FALSE)
    }
  })
  vc <- rvinecopulib::vinecop_dist(pcs, rvinecopulib::dvine_structure(seq_len(n)))
  log(rvinecopulib::dvinecop(matrix(u, 1), vc))
}

pair <- function(family, rotation = 0, parameters) {
  list(family = family, rotation = rotation, parameters = parameters)
}
