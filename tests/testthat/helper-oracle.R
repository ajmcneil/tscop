# Oracle tests compare against the legacy tscopula package. It is only a
# numerical reference: no tscopula code is used in tscop.
skip_if_no_oracle <- function() {
  testthat::skip_if_not_installed("tscopula")
}
