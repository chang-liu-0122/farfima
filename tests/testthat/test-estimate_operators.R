test_that("estimate_operators gives the operators of the fit and draws them", {
  fit = farfima(X[, 1:100], d = 0.3, order = c(1, 1), trim = 20, verbose = FALSE)
  Y = frac_diff(X[, 1:100] - rowMeans(X[, 1:100]), d = 0.3)[, -(1:20)]
  op = estimate_operators(Y, order = c(1, 1))
  expect_s3_class(fit$operators, "farfima_operators")
  expect_equal(op$phi, fit$phi)
  expect_equal(op$theta, fit$theta)
  expect_equal(op$residuals, fit$residuals, ignore_attr = TRUE)
  expect_equal(dim(op$phi), c(21, 21, 1))
  expect_equal(op$fitted[, -1] + op$residuals[, -1], Y[, -1], ignore_attr = TRUE)
  expect_output(print(op), "Norms: +AR lag 1 = .*, MA lag 1 = ")
  o = select_order(Y, p = 0:1, q = 0:1)
  expect_equal(estimate_operators(Y, order = o)$order, c(p = o$p, q = o$q))
  expect_error(estimate_operators(Y, order = c(0, 0)), "p \\+ q > 0")
  grDevices::pdf(NULL)
  expect_silent(plot(op))
  expect_silent(plot(estimate_operators(Y, order = c(3, 2)), main = "FARMA(3, 2)"))
  grDevices::dev.off()
})

test_that("estimate_operators warns when the fitted operators are unstable", {
  local_mocked_bindings(farma_fit_diagnostics = function(fit, Y_train) {
    list(ar_radius = 0.5, ma_radius = 1.2, residual_ratio = 3)
  })
  expect_warning(estimate_operators(X, order = c(1, 1)), "FARMA\\(1, 1\\) is unstable; try smaller p and q")
})
