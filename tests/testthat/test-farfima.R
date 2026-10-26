test_that("farfima returns the model and the calibration table", {
  fit = farfima(X[, 1:100], X[, 101:130], d = 0.3, order = c(1, 0), h = 3, trim = 20, verbose = FALSE)
  expect_named(fit, c("d", "order", "phi", "theta", "fitted", "fitted_differenced", "residuals",
                      "memory", "order_selection", "operators", "calibration"))
  expect_named(fit$calibration, c("level", "h", "SR_zeta", "NS_quantile", "AC_eta", "n_origin"))
  expect_identical(fit$order, c(p = 1L, q = 0L))
})

test_that("an unstable final fit gives a warning and is shown when printed", {
  local_mocked_bindings(farma_fit_diagnostics = function(fit, Y_train) {
    list(ar_radius = 0.5, ma_radius = 1.2, residual_ratio = 3)
  })
  expect_warning({
    fit = farfima(X[, 1:100], d = 0.3, order = c(1, 1), trim = 20, verbose = FALSE)
  }, "FARFIMA\\(1, 0.3, 1\\) is unstable; try smaller p and q")
  expect_output(print(fit), "Stability: unstable \\(MA part not invertible\\)")
  expect_equal(attr(fit, "stability")$ma_radius, 1.2)
})

test_that("failed validation origins are skipped with a warning", {
  real = forecast_at_origin
  local_mocked_bindings(forecast_at_origin = function(X, t0, ...) {
    if (t0 == 105) stop(unstable(1.1))
    real(X, t0, ...)
  })
  expect_warning({
    fit = farfima(X[, 1:100], X[, 101:130], d = 0.3, order = c(1, 0), h = 3, trim = 20, verbose = FALSE)
  }, "1 of 30 validation forecasts failed and were skipped \\(origins 105\\)")
  expect_equal(fit$calibration$n_origin[fit$calibration$level == 0.8], c(29, 28, 27))
  expect_output(print(fit), "1 of 30 validation forecasts failed and were skipped")
})

test_that("the bands cannot be calibrated when every validation origin fails", {
  local_mocked_bindings(forecast_at_origin = function(X, t0, ...) stop(unstable(1.5)))
  expect_error(suppressWarnings(
    farfima(X[, 1:100], X[, 101:130], d = 0.3, order = c(1, 0), h = 3, trim = 20, verbose = FALSE)
  ), "No validation forecast succeeded")
})

test_that("the PID and AcMCP learning rates are calibrated on the validation sample", {
  for (sym in c(TRUE, FALSE)) {
    fit = farfima(X[, 1:100], X[, 101:130], d = 0.3, order = c(1, 0), h = 3, trim = 20,
                   band = c("NS", "AC", "PID", "AcMCP"), symmetric = sym, verbose = FALSE)
    expect_named(fit$calibration, c("level", "h", "SR_zeta", "NS_quantile", "AC_eta", "PID_eta", "AcMCP_eta", "n_origin"))
    expect_true(all(fit$calibration$PID_eta %in% c(0, 0.005, 0.01, 0.02, 0.05, 0.1)))
  }
})

test_that("the scorecaster is shown when printed and must be naive or theta", {
  fit = farfima(X[, 1:100], X[, 101:130], d = 0.3, order = c(1, 0), h = 3, trim = 20, band = c("NS", "PID"),
                 scorecast = "theta", verbose = FALSE)
  expect_output(print(fit), "PID with theta scorecaster")
  expect_error(farfima(X[, 1:100], X[, 101:130], d = 0.3, order = c(1, 0), h = 3, trim = 20, scorecast = "arima"), "should be one of")
})
