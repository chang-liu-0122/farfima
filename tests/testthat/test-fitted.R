test_that("fitted and residuals give back the training curves", {
  fit = farfima(X[, 1:100], d = 0.3, order = c(1, 0), trim = 20, verbose = FALSE)
  expect_equal(dim(residuals(fit)), c(21, 80))
  expect_true(all(is.na(residuals(fit)[, 1])))
  expect_equal((fitted(fit) + residuals(fit))[, -1], X[, 21:100][, -1], ignore_attr = TRUE)
  expect_equal(dim(fitted(fit, scale = "differenced")), c(21, 80))
  bench = farfima(X[, 1:100], d = 1, method = "far", trim = 20, verbose = FALSE)
  expect_error(fitted(bench), "only available for the FARMA method")
})
