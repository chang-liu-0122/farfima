test_that("frac_diff with -d undoes frac_diff with d", {
  y = frac_diff(X, 0.3)
  expect_equal(dim(y), dim(X))
  expect_equal(frac_diff(y, -0.3), X, ignore_attr = TRUE)
  expect_equal(frac_diff(X, 0), X, ignore_attr = TRUE)
  expect_equal(frac_diff(X, 1)[, -1], X[, -1] - X[, -ncol(X)], ignore_attr = TRUE)
})
