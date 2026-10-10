# Automatic selection of non-Gaussian pair copulas for d-vine models.

.auto_choices <- c("gumbel", "clayton", "frank", "joe")
.auto_allowed <- c(.auto_choices, "ast")

# Information criterion of a fit.
dv_ic <- function(loglik, npar, nobs, criterion) {
  switch(criterion,
    AIC = -2 * loglik + 2 * npar,
    BIC = -2 * loglik + log(nobs) * npar,
    AICc = -2 * loglik + 2 * npar + 2 * npar * (npar + 1) / (nobs - npar - 1)
  )
}

dv_fit_ic <- function(fit, criterion) dv_ic(fit@loglik, fit@npar, fit@nobs, criterion)

# Free vector of a fully specified tied model (ARMA coefficients, auxiliary
# parameters).
dv_model_z <- function(model) {
  lay <- dv_layout(model)
  aux <- vapply(model@pairs[lay$aux_lags], function(sp) sp@aux, numeric(1))
  c(lay$ops$to_free(model@kpacf), dv_aux_to_free(aux, lay$aux_fam),
    if (lay$has_ivt) dv_vt_to_free(model@vt, lay$vtpar))
}

# Candidate pair copula specifications at a lag where Kendall's tau is
# positive or negative: only the rotations for that sign are tried (0 and 180
# for positive dependence, 90 and 270 for negative); the other sign keeps its
# default rotation.
dv_candidates <- function(choices, aux, positive) {
  out <- list()
  for (fam in c(choices, aux)) {
    if (fam %in% .pc_radial) {
      out[[length(out) + 1L]] <- pc(fam)
    } else {
      for (r in if (positive) c(0, 180) else c(90, 270)) {
        out[[length(out) + 1L]] <- if (positive) pc(fam, pos = r, neg = 90) else pc(fam, pos = 0, neg = r)
      }
    }
  }
  out
}

# Score candidates at lag k by a single log-likelihood evaluation at the
# current parameters (the candidate replaces the pair copula at lag k, all
# other pair copulas unchanged); for t and BB1 the auxiliary parameter is found
# by a one-dimensional search. Returns a data frame with the information
# criterion, log-likelihood and auxiliary value for each candidate.
dv_screen <- function(u, best_model, k, cands, criterion, npar_best) {
  lay <- dv_layout(best_model)
  K_use <- max(best_model@K, k)
  pairs_base <- dv_make_pairs(dv_model_z(best_model), lay, K_use)
  tk <- gc_tau(best_model@kpacf, K_use)[k]
  n <- length(u)
  ll_at <- function(sp, a) {
    pr <- dv_resolve_vals(sp@family, sp@pos, sp@neg, a, tk)
    if (is.null(pr)) return(-Inf)
    pb <- pairs_base
    pb[[k]] <- pr
    ll <- tryCatch(dv_loglik(u, pb), error = function(err) -Inf)
    if (is.finite(ll)) ll else -Inf
  }
  rows <- lapply(cands, function(sp) {
    fam <- sp@family
    a <- NA_real_
    ll <- -Inf
    if (fam == "t") {
      o <- stats::optimize(function(x) { v <- ll_at(sp, 2 + exp(x)); if (is.finite(v)) -v else 1e10 },
                           c(log(0.05), log(200)))
      a <- 2 + exp(o$minimum)
      ll <- -o$objective
    } else if (fam == "bb1") {
      dmax <- 1 / (1 - min(abs(tk), 0.999))
      if (dmax > 1 + 1e-4) {
        o <- stats::optimize(function(d) { v <- ll_at(sp, d); if (is.finite(v)) -v else 1e10 },
                             c(1, dmax * (1 - 1e-6)))
        a <- o$minimum
        ll <- -o$objective
      }
    } else {
      ll <- ll_at(sp, NA_real_)
    }
    if (ll <= -1e9) ll <- -Inf
    p <- npar_best + as.integer(fam %in% .pc_aux_families)
    data.frame(ic = if (is.finite(ll)) dv_ic(ll, p, n, criterion) else Inf, loglik = ll, aux = a)
  })
  do.call(rbind, rows)
}

# Sign of Kendall's tau at lag k of a model's ARMA part (TRUE: non-negative).
dv_tau_positive <- function(model, k) {
  gc_tau(model@kpacf, k)[k] >= 0
}

# Starting model for the search: a fitted d-vine of the same data.
dv_auto_start <- function(U, start, tautol, maxlag, pseudo, control) {
  if (is(start, "dvinecopfit")) return(tsfit(U, start@model, pseudo = pseudo, control = control))
  if (is(start, "dvinecop")) return(tsfit(U, start, pseudo = pseudo, control = control))
  if (is(start, "armacopfit") || is(start, "arfactcopfit")) start <- start@model
  if (is(start, "gausscop")) {
    return(tsfit(U, dvinecop(start, tautol = tautol, maxlag = maxlag), pseudo = pseudo,
                 control = control))
  }
  stop("'start' must be an armacop (or arfactcop) model or fit, or a dvinecop or dvinecopfit object")
}

#' Automatic selection of non-Gaussian pair copulas
#'
#' Starting from a Gaussian ARMA copula (or an existing d-vine fit),
#' systematically replaces Gaussian pair copulas by non-Gaussian ones, lag by
#' lag, keeping the substitution that most improves an information criterion.
#' It is the tool for users with no prior idea of which copulas to use;
#' [tsfit()] fits a fully specified model.
#'
#' The search is greedy and sequential. At lag `k` the sign of the Kendall
#' partial autocorrelation at the current estimates decides which rotations
#' are tried (0 and 180 degrees for positive dependence, 90 and 270 for
#' negative), so the cost is linear in the number of lags. With
#' `screen = TRUE` (the default) each candidate family and rotation is scored
#' by a single log-likelihood evaluation at the current ARMA coefficients (for
#' `t` and BB1 copulas after a one-dimensional search over the auxiliary
#' parameter); only the best `ntop` candidates are refitted jointly with
#' [tsfit()]. With `screen = FALSE` every candidate is refitted. If the sign of
#' Kendall's tau at the lag changes in a refit, the rotations for the new sign
#' are also tried. A substitution is accepted if it lowers the criterion by more than
#' `accepttol`; the lag stays Gaussian otherwise. Lags at which Kendall's tau
#' is below the effective-lag tolerance `tautol` are skipped (every family is
#' then the independence copula). The search stops after `nstrike` consecutive
#' lags whose improvement is below `ICtol`, or when lag `nreplace` has been
#' examined. Standard errors are not computed; run [tsfit()] with `se` on the
#' selected model (which ignores the uncertainty from the selection itself).
#'
#' @param U complete copula data in (0, 1).
#' @param start the starting model: an [armacop-class] model or fit (the
#'   Gaussian d-vine with this ARMA kpacf is fitted first), or a tied
#'   [dvinecop-class] model or fit (to continue a search).
#' @param nreplace the largest lag at which a substitution is examined (the
#'   search examines lags after any explicit pair copulas already present, up
#'   to `nreplace`).
#' @param choices one-parameter families to try: any of `"gumbel"`,
#'   `"clayton"`, `"frank"`, `"joe"` (the default) and `"ast"` (the absolute
#'   spherical t copula, not in the default because it is about ten times
#'   slower to evaluate).
#' @param aux two-parameter families to add to the candidates: `"t"` and/or
#'   `"bb1"`; `NULL` for none.
#' @param criterion `"AIC"`, `"BIC"` or `"AICc"`; the number of parameters
#'   counts the ARMA coefficients and every auxiliary parameter.
#' @param ICtol,nstrike stop after `nstrike` consecutive lags whose
#'   improvement in the criterion is less than `ICtol`.
#' @param screen logical: screen candidates at fixed parameters before refitting.
#' @param ntop number of screened candidates refitted jointly.
#' @param accepttol the smallest improvement in the criterion for which a
#'   substitution is accepted (default 0: any improvement, as in the original
#'   `auto_dvine`). With many candidate families some chance improvements are
#'   to be expected on Gaussian data; a positive value such as `ICtol` guards
#'   against them.
#' @param pseudo `TRUE` if `U` are pseudo-observations; passed to [tsfit()].
#' @param tautol,maxlag effective lag settings for new models (ignored if
#'   `start` is already a d-vine).
#' @param control list of control parameters for [stats::optim()].
#' @param verbose logical: print a line per lag.
#'
#' @return A [dvinecopfit-class] object for the selected model, with the
#'   selection path in its `path` slot.
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' truth <- dvinecop(armacop(ar = c(0.6, -0.2)), pairs = list(pc("gumbel"), pc("frank")))
#' u <- tssim(truth, 800)
#' gauss <- tsfit(u, armacop(order = c(2, 0)), pseudo = FALSE)
#' autodvinecop(u, gauss, nreplace = 4, pseudo = FALSE)
#' }
autodvinecop <- function(U, start, nreplace = 10, choices = .auto_choices, aux = NULL,
                         criterion = "AIC", ICtol = 0.5, nstrike = 3, screen = TRUE,
                         ntop = 2, accepttol = 0, pseudo = TRUE, tautol = 1e-4, maxlag = Inf,
                         control = list(), verbose = TRUE) {
  criterion <- match.arg(criterion, c("AIC", "BIC", "AICc"))
  choices <- tolower(choices)
  if (!all(choices %in% .auto_allowed)) {
    stop("'choices' must be among ", paste(.auto_allowed, collapse = ", "))
  }
  if (!is.null(aux)) {
    aux <- tolower(aux)
    if (!all(aux %in% .pc_aux_families)) stop("'aux' must be among ", paste(.pc_aux_families, collapse = ", "))
  }
  if (length(choices) + length(aux) == 0L) stop("no candidate families")
  if (!is.numeric(nreplace) || length(nreplace) != 1L || nreplace < 1 || nreplace != round(nreplace)) {
    stop("'nreplace' must be a positive whole number")
  }
  if (!is.numeric(nstrike) || nstrike < 1) stop("'nstrike' must be at least 1")
  u <- check_complete_copula_data(U, min_n = 20L)
  n <- length(u)
  if (is(start, "dvinecop") && start@mode != "tied") stop("'start' must be a tied d-vine model")
  if (is(start, "dvinecopfit") && start@model@mode != "tied") stop("'start' must be a tied d-vine model")

  best <- dv_auto_start(U, start, tautol, maxlag, pseudo, control)
  K0 <- length(best@model@pairs)
  pad <- 0L
  strike <- 0L
  path <- list()
  say <- function(...) if (verbose) cat(..., "\n", sep = "")
  say("Starting model: ", criterion, " = ", format(dv_fit_ic(best, criterion), digits = 7),
      ", ", K0, " explicit pair copula(s), effective lag ", best@model@K)

  ks <- if (K0 < nreplace) (K0 + 1L):as.integer(nreplace) else integer(0)
  if (length(ks) == 0L) say("Nothing to do: the model already has ", K0, " explicit pair copulas")
  for (k in ks) {
    cur <- best@model
    ic0 <- dv_fit_ic(best, criterion)
    positive <- dv_tau_positive(cur, k)
    tau_k <- {
      gc_tau(cur@kpacf, k)[k]
    }
    if (abs(tau_k) <= cur@tautol) {
      # beyond the kpacf's effective range every family is the independence
      # copula: nothing can be gained, and the lag stays Gaussian
      pad <- pad + 1L
      path[[length(path) + 1L]] <- data.frame(
        lag = k, choice = "gauss", tau = tau_k, screened = NA_real_, criterion = ic0,
        improvement = 0, accepted = FALSE, nfit = 0L, stringsAsFactors = FALSE
      )
      say("lag ", k, ": |tau| = ", format(abs(tau_k), digits = 3), " is below tautol; skipped")
      strike <- strike + 1L
      if (strike >= nstrike) {
        say("Stopped: ", nstrike, " consecutive lags with improvement below ", ICtol)
        break
      }
      next
    }
    cands <- dv_candidates(choices, aux, positive)
    screened_ic <- NA_real_
    if (screen) {
      sc <- dv_screen(u, cur, k, cands, criterion, best@npar)
      ord <- order(sc$ic)
      ord <- ord[is.finite(sc$ic[ord])]
      ord <- utils::head(ord, ntop)
      screened_ic <- if (length(ord) > 0L) sc$ic[ord[1L]] else NA_real_
      sel <- lapply(ord, function(i) {
        sp <- cands[[i]]
        pc(sp@family, sp@pos, sp@neg, if (is.na(sc$aux[i])) NULL else sc$aux[i])
      })
    } else {
      sel <- cands
    }
    base_pairs <- c(cur@pairs, rep(list(pc("gauss")), pad))
    fit_with <- function(sp) {
      m <- dvinecop(cur@kpacf, pairs = c(base_pairs, list(sp)), tail = cur@tail,
                    tautol = cur@tautol, maxlag = cur@maxlag,
                    vt = if (length(cur@vt)) cur@vt else NULL, vtpar = cur@vtpar)
      tryCatch(suppressWarnings(tsfit(u, m, pseudo = pseudo, control = control)),
               error = function(err) NULL)
    }
    fits <- lapply(sel, fit_with)
    tried <- sel
    # if the sign of tau at the lag flipped in a refit, try the other sign's rotations
    extra <- list()
    for (i in seq_along(fits)) {
      f <- fits[[i]]
      sp <- sel[[i]]
      if (is.null(f) || sp@family %in% .pc_radial) next
      pos_now <- dv_tau_positive(f@model, k)
      if (pos_now != positive) {
        for (r in if (pos_now) c(0, 180) else c(90, 270)) {
          v <- if (pos_now) pc(sp@family, r, sp@neg, if (is.na(sp@aux)) NULL else sp@aux)
               else pc(sp@family, sp@pos, r, if (is.na(sp@aux)) NULL else sp@aux)
          if (!any(vapply(tried, function(t) identical(t, v), logical(1)))) extra[[length(extra) + 1L]] <- v
        }
      }
    }
    if (length(extra) > 0L) {
      fits <- c(fits, lapply(extra, fit_with))
      tried <- c(tried, extra)
    }
    ok <- !vapply(fits, is.null, logical(1))
    accepted <- FALSE
    ic_new <- ic0
    choice <- "gauss"
    nfit <- sum(ok)
    if (any(ok)) {
      ics <- vapply(fits[ok], dv_fit_ic, numeric(1), criterion = criterion)
      j <- which.min(ics)
      if (ics[j] < ic0 - accepttol) {
        best <- fits[ok][[j]]
        ic_new <- ics[j]
        accepted <- TRUE
        pad <- 0L
        choice <- dv_explicit_labels(best@model)[k]
      }
    }
    if (!accepted) pad <- pad + 1L
    improvement <- ic0 - ic_new
    path[[length(path) + 1L]] <- data.frame(
      lag = k, choice = choice, tau = tau_k, screened = screened_ic, criterion = ic_new,
      improvement = improvement, accepted = accepted, nfit = nfit, stringsAsFactors = FALSE
    )
    say("lag ", k, ": tau = ", format(tau_k, digits = 3), "; ",
        if (accepted) paste0("selected ", choice) else "no improvement (stays Gaussian)",
        "; ", criterion, " = ", format(ic_new, digits = 7), " (improvement ",
        format(improvement, digits = 3), ")")
    if (improvement < ICtol) strike <- strike + 1L else strike <- 0L
    if (strike >= nstrike) {
      say("Stopped: ", nstrike, " consecutive lags with improvement below ", ICtol)
      break
    }
  }
  best@path <- if (length(path) > 0L) do.call(rbind, path) else data.frame()
  best
}
