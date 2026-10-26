test_that("every memory estimator returns a finite value", {
  for (m in c("lw", "dfa", "lpwn1", "lpwn2")) {
    expect_true(is.finite(estimate_memory(X, method = m)$d))
  }
  bc = estimate_memory(X, method = "bc_lpwn2", B = 3, burn_in = 20, seed = 2)
  expect_named(bc, c("d", "scores", "d_initial", "bias", "bootstrap"))
  expect_named(estimate_memory(X, method = "lw"), c("d", "scores"))
  expect_equal(bc$d, bc$d_initial - bc$bias)
  expect_equal(bc$d, estimate_memory(X, method = "bc_lpwn2", B = 3, burn_in = 20, seed = 2)$d)
})
