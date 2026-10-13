test_that("Joe and Frank tau-parameter relations agree with rvinecopulib and are smooth", {
  th <- c(1.0001, 1.2, 1.9999, 2.0001, 3, 7, 20, 30)
  expect_equal(joe_tau(th), vapply(th, function(x) as.numeric(rvinecopulib::par_to_ktau("joe", 0, x)), 1), tolerance = 1e-10)
  expect_equal(joe_tau(2), 2 - pi^2 / 6, tolerance = 1e-14)        # removable singularity at theta = 2
  th <- c(0.001, 0.1, 1, 1.99, 2.01, 3, 5, 10, 30, 35)
  expect_equal(frank_tau(th), vapply(th, function(x) as.numeric(rvinecopulib::par_to_ktau("frank", 0, x)), 1), tolerance = 1e-10)
  expect_equal(frank_tau(-th), -frank_tau(th))
  h <- 1e-6
  for (x in c(1.3, 2, 4, 15)) {
    expect_equal(joe_dtau(x), (joe_tau(x + h) - joe_tau(x - h)) / (2 * h), tolerance = 1e-7)
  }
  for (x in c(0.5, 1.99, 2.01, 4, 15)) {
    expect_equal(frank_dtau(x), (frank_tau(x + h) - frank_tau(x - h)) / (2 * h), tolerance = 1e-7)
  }
})

test_that("the inverses are accurate, vectorised, and the compiled versions equal the R ones", {
  tau <- c(1e-9, 1e-6, 1e-3, 0.05, 0.2, 0.5, 0.8, 0.9, 0.99)
  pj <- joe_par(tau)
  expect_equal(joe_tau(pj), tau, tolerance = 1e-9)
  expect_equal(joe_par_fast(tau), pj, tolerance = 1e-9)
  expect_equal(joe_par(0), 1)
  tf <- c(-tau, tau)
  pf <- frank_par(tf)
  expect_equal(frank_tau(pf), tf, tolerance = 1e-12)
  expect_equal(frank_par_fast(tf), pf, tolerance = 1e-12)
  expect_equal(frank_par(0), 1e-10)
  expect_equal(frank_par_fast(c(0, -1e-15)), c(1e-10, -1e-10))
  a <- joe_tau_d(c(1.5, 2, 9)); b <- frank_tau_d(c(-3, 0.5, 2.5))
  expect_equal(a$tau, joe_tau(c(1.5, 2, 9)), tolerance = 1e-14)
  expect_equal(a$dtau, joe_dtau(c(1.5, 2, 9)), tolerance = 1e-12)
  expect_equal(b$tau, frank_tau(c(-3, 0.5, 2.5)), tolerance = 1e-14)
  expect_equal(b$dtau, frank_dtau(c(-3, 0.5, 2.5)), tolerance = 1e-12)
  # the numerical inverse of rvinecopulib is only accurate to about 1e-7 in the parameter
  old <- vapply(tau[2:6], function(t) as.numeric(rvinecopulib::ktau_to_par("frank", t)), 1)
  expect_equal(frank_par_fast(tau[2:6]), old, tolerance = 1e-5)
})

test_that("dv_ktau_to_par uses the accurate inverses", {
  tau <- c(0.03, 0.3, 0.7)
  expect_equal(rvinecopulib::par_to_ktau("joe", 0, dv_ktau_to_par("joe", tau[2])), tau[2], tolerance = 1e-12)
  expect_equal(vapply(dv_ktau_to_par("frank", tau), function(p) as.numeric(rvinecopulib::par_to_ktau("frank", 0, p)), 1),
               tau, tolerance = 1e-12)
})
