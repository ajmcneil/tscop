# Models with BB1 pair copulas (feasible auxiliary parameters: delta < 1 / (1 - |tau|) at the lag)
bb1_models <- function() {
  kp <- arfactcop(lambda = 0.9, w = 0.5)
  list(
    lag1 = sdvinecop(kp, pairs = list(pc("bb1", aux = 1.2)), tail = pc("joe"), maxlag = 10, tautol = 1e-4),
    two = sdvinecop(kp, pairs = list(pc("bb1", aux = 1.25), pc("bb1", aux = 1.1), pc("gumbel")), tail = pc("clayton"),
                    maxlag = 10, tautol = 1e-4),
    rot = sdvinecop(armacop(ar = c(0.6, -0.3), ma = 0.2),
                    pairs = list(pc("bb1", aux = 1.4), pc("bb1", aux = 1.2, pos = 180, neg = 270)), tail = pc("frank"),
                    maxlag = 10, tautol = 1e-4)
  )
}
