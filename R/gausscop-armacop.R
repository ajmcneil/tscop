# The gausscop interface for armacop.

setMethod("gc_is_specified", "armacop", function(G) arma_is_specified(G))

setMethod("gc_describe", "armacop", function(G) {
  paste0(G@name, if (!is.na(G@period)) paste0(" with period ", G@period) else "")
})

setMethod("gc_pacf", "armacop", function(G, K) {
  e <- arma_coef_expanded(G)
  arma_pacf(e$ar, e$ma, K)
})

setMethod("gc_effective_lag", "armacop", function(G, K0, tautol, maxlag) {
  e <- arma_coef_expanded(G)
  dv_effective_lag(e$ar, e$ma, K0, tautol, maxlag, G@period)
})

setMethod("gc_structure", "armacop", function(G) list(class = "armacop", order = G@order, period = G@period))

setMethod("gc_prefit", "armacop", function(G, U, control) armacop_mle(U, G, NULL, control)$model)

setMethod("gc_from_coef", "armacop", function(G, v) {
  pars <- split_blocks(as.numeric(v), G@order)
  armacop(ar = pars$ar, ma = pars$ma, sar = pars$sar, sma = pars$sma, period = G@period)
})

setMethod("gc_free_ops", "armacop", function(G) {
  ord <- G@order
  period <- G@period
  list(
    nfree = sum(ord),
    zbound = 6,
    to_free = function(G2) pars_to_free(G2@pars),
    pacf = function(z, K) {
      e <- free_to_expanded(z, ord, period)
      arma_pacf(e$ar, e$ma, K)
    },
    effective_lag = function(z, K0, tautol, maxlag) {
      e <- free_to_expanded(z, ord, period)
      dv_effective_lag(e$ar, e$ma, K0, tautol, maxlag, period)
    },
    coef = function(z) free_to_coef(z, ord),
    build = function(z) {
      b <- split_blocks(z, ord)
      armacop(ar = free_to_ar(b$ar), ma = free_to_ma(b$ma), sar = free_to_ar(b$sar),
              sma = free_to_ma(b$sma), period = period)
    }
  )
})
