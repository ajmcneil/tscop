# Models with absolute spherical t (ast) pair copulas: plain, as an IVT base, and mixed with other families
ast_models <- function() {
  kp <- arfactcop(lambda = 0.95, w = 0.5)
  list(
    plain = sdvinecop(kp, tail = pc("ast"), maxlag = 12, tautol = 1e-6),
    ivt = sdvinecop(kp, tail = pc("ivt", base = "ast"), vt = c(0.45, 0.55), maxlag = 12, tautol = 1e-6),
    ivt_asym = sdvinecop(kp, tail = pc("ivt", base = "ast"), vt = c(0.3, 0.7), maxlag = 12, tautol = 1e-6),
    ast_then_gumbel = sdvinecop(kp, pairs = rep(pc("ast"), 3), tail = pc("gumbel"), maxlag = 12, tautol = 1e-6),
    gumbel_then_ast = sdvinecop(kp, pairs = list(pc("gumbel")), tail = pc("ast"), maxlag = 12, tautol = 1e-6),
    ivt_ast_then_joe = sdvinecop(kp, pairs = rep(pc("ivt", base = "ast"), 3), tail = pc("ivt", base = "joe"), vt = c(0.45, 0.55),
                                 maxlag = 12, tautol = 1e-6),
    joe_then_ivt_ast = sdvinecop(kp, pairs = rep(pc("ivt", base = "joe"), 2), tail = pc("ivt", base = "ast"), vt = c(0.45, 0.55),
                                 maxlag = 12, tautol = 1e-6),
    mixed = sdvinecop(kp, pairs = list(pc("gauss"), pc("ast"), pc("gauss"), pc("ast"), pc("ivt", base = "ast")),
                      tail = pc("ivt", base = "ast"), vt = c(0.4, 0.6), maxlag = 12, tautol = 1e-6),
    adjacent = sdvinecop(kp, pairs = c(rep(list(pc("ast")), 2), rep(list(pc("ivt", base = "ast")), 2)), tail = pc("ast"),
                         vt = c(0.4, 0.6), maxlag = 12, tautol = 1e-6)
  )
}

# A model the kernels do not support: ast pairs with a negative tau at lag 2 are rotated (rotation 90)
unsupported_model <- function() {
  m <- sdvinecop(armacop(ar = c(0.5, -0.3)), pairs = list(pc("ast"), pc("ast")), tail = pc("gauss"), maxlag = 6)
  m@K <- 6L
  m
}
