test_that("plot draws the fitted and residual curves", {
  fit = farfima(X[, 1:100], d = 0.3, order = c(1, 0), trim = 20, verbose = FALSE)
  grDevices::pdf(NULL)
  expect_silent(plot(fit))
  expect_silent(plot(fit, type = "residuals", main = "Residuals"))
  grDevices::dev.off()
  bench = farfima(X[, 1:100], d = 1, method = "far", trim = 20, verbose = FALSE)
  expect_error(plot(bench), "only available for the FARMA method")
})
