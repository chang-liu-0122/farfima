test_that("plot_curves draws the rainbow plot and the periods", {
  Y = X[, 1:40]
  colnames(Y) = 1961:2000
  grDevices::pdf(NULL)
  expect_silent(plot_curves(Y))
  expect_silent(plot_curves(Y, type = "periods", periods = 4))
  expect_silent(plot_curves(Y, type = "periods", periods = c(1961, 1975, 1990), main = "Simulated curves"))
  expect_silent(plot_curves(X[, 1:30], argvals = seq(0, 1, length.out = 21)))
  expect_identical(plot_curves(Y, type = "periods"), Y)
  grDevices::dev.off()
  expect_error(plot_curves(Y, type = "surface"), "should be one of")
})
