# stats::ARMAacf rejects the empty model; white noise has acf (1, 0, 0, ...).
acf_ref <- function(ar, ma, lag.max) {
  if (length(ar) == 0 && length(ma) == 0) return(c(1, numeric(lag.max)))
  as.numeric(ARMAacf(ar = ar, ma = ma, lag.max = lag.max))
}
