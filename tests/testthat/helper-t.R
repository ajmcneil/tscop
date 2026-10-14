# Models with t copula pairs (the degrees of freedom are the auxiliary parameter), also mixed with ast and others
t_models <- function() {
  kp <- arfactcop(lambda = 0.9, w = 0.5)
  list(
    lag1 = sdvinecop(kp, pairs = list(pc("t", aux = 6)), tail = pc("joe"), maxlag = 10, tautol = 1e-4),
    two = sdvinecop(kp, pairs = list(pc("t", aux = 6), pc("t", aux = 9), pc("gumbel")), tail = pc("clayton"), maxlag = 10, tautol = 1e-4),
    run_then_ast = sdvinecop(kp, pairs = list(pc("t", aux = 5), pc("t", aux = 7), pc("t", aux = 12)), tail = pc("ast"), maxlag = 10, tautol = 1e-4),
    negative = sdvinecop(armacop(ar = c(0.6, -0.3), ma = 0.2), pairs = list(pc("t", aux = 6), pc("t", aux = 8), pc("t", aux = 4)),
                         tail = pc("frank"), maxlag = 10, tautol = 1e-4),
    after_gauss = sdvinecop(kp, pairs = list(pc("gauss"), pc("t", aux = 6), pc("t", aux = 5), pc("ivt", base = "joe")),
                            tail = pc("ivt", base = "joe"), vt = c(0.45, 0.55), maxlag = 10, tautol = 1e-4)
  )
}
