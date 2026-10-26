test_that("predict returns bands and an optional test evaluation", {
  fit = farfima(X[, 1:100], X[, 101:130], d = 0.3, order = c(1, 0), h = 3, trim = 20, verbose = FALSE)
  fc = predict(fit, verbose = FALSE)
  expect_equal(dim(fc$point_forecast), c(21, 3))
  expect_null(fc$evaluation)
  b = fc$bands$SR[["80%"]]
  expect_true(all(b$lower <= fc$point_forecast & fc$point_forecast <= b$upper))
  fc2 = predict(fit, test = X[, 131:160], verbose = FALSE)
  expect_named(fc2, c("point_forecast", "bands", "evaluation"))
  expect_named(fc2$evaluation, c("point", "interval", "wis"))
  expect_equal(nrow(fc2$evaluation$point), 3)
  expect_setequal(unique(fc2$evaluation$interval$band), c("SR", "NS", "AC"))
  expect_named(fc2$evaluation$interval, c("band", "level", "h", "ECP", "CPD", "MIS"))
  r = Filter(function(z) length(z$target_test_curve) >= 2 && z$target_test_curve[2] == 5, attr(fc2, "test_forecasts"))[[1]]
  expect_equal(r$observed[, 2], X[, 130 + 5])
  expect_error(predict(fit, h = 4, verbose = FALSE), "calibrated for h <= 3")
  expect_error(predict(fit, level = 0.9), "set in farfima")
})

test_that("predict keeps only the bands chosen in farfima", {
  fit_ac = farfima(X[, 1:100], X[, 101:130], d = 0.3, order = c(1, 0), h = 3, band = "AC", trim = 20, verbose = FALSE)
  expect_equal(names(predict(fit_ac, verbose = FALSE)$bands), "AC")
  expect_error(predict(fit_ac, band = "SR", verbose = FALSE), "not chosen in farfima")
})

test_that("the FAR(1) benchmark forecasts the test sample", {
  bench = farfima(X[, 1:100], X[, 101:130], d = 1, method = "far", h = 3, trim = 20, verbose = FALSE)
  fb = predict(bench, test = X[, 131:160], verbose = FALSE)
  expect_true(all(is.finite(fb$evaluation$point$RMSPE)))
})

test_that("samples out of time order give a warning", {
  Y = X
  colnames(Y) = 1850 + seq_len(ncol(Y))
  fit = farfima(Y[, 1:100], d = 0.3, order = c(1, 0), trim = 20, verbose = FALSE)
  expect_warning(predict(fit, test = Y[, 131:160], h = 3, verbose = FALSE), "must follow each other")
  expect_silent(predict(fit, test = Y[, 101:130], h = 3, verbose = FALSE))
})

test_that("a fit without validation gives point forecasts only", {
  fit = farfima(X[, 1:100], d = 0.3, order = c(1, 0), trim = 20, verbose = FALSE)
  fc = predict(fit, test = X[, 101:130], h = 3, verbose = FALSE)
  expect_length(fc$bands, 0)
  expect_null(fc$evaluation$interval)
  expect_equal(nrow(fc$evaluation$point), 3)
})

test_that("failed test origins are skipped with a warning", {
  fit = farfima(X[, 1:100], X[, 101:130], d = 0.3, order = c(1, 0), h = 3, trim = 20, verbose = FALSE)
  real = forecast_at_origin
  local_mocked_bindings(forecast_at_origin = function(X, t0, ...) {
    if (t0 == 140) stop(unstable(1.1))
    real(X, t0, ...)
  })
  expect_warning({
    fc = predict(fit, test = X[, 131:160], verbose = FALSE)
  }, "1 of 30 test forecasts failed and were skipped \\(origins 140\\)")
  expect_equal(names(attr(fc, "failures")), "140")
  expect_true(all(is.finite(fc$evaluation$point$RMSPE)))
})

test_that("an unstable refit for the final forecast gives a warning but still forecasts", {
  fit = farfima(X[, 1:100], X[, 101:130], d = 0.3, order = c(1, 0), h = 3, trim = 20, verbose = FALSE)
  real = forecast_at_origin
  local_mocked_bindings(forecast_at_origin = function(X, t0, ..., check = TRUE) {
    if (t0 == ncol(X) && check) stop(unstable(1.05))
    real(X, t0, ..., check = check)
  })
  expect_warning({
    fc = predict(fit, verbose = FALSE)
  }, "forecast beyond the data comes from an unstable FARMA fit")
  expect_true(all(is.finite(fc$point_forecast)))
})

test_that("the PID and AcMCP bands run on the test sample", {
  for (sym in c(TRUE, FALSE)) {
    fit = farfima(X[, 1:100], X[, 101:130], d = 0.3, order = c(1, 0), h = 3, trim = 20,
                   band = c("NS", "AC", "PID", "AcMCP"), symmetric = sym, verbose = FALSE)
    fc = predict(fit, test = X[, 131:160], verbose = FALSE)
    expect_setequal(unique(fc$evaluation$interval$band), c("NS", "AC", "PID", "AcMCP"))
    expect_true(all(is.finite(fc$evaluation$interval$MIS)))
    for (m in c("AC", "PID", "AcMCP")) {
      b = fc$bands[[m]][["95%"]]
      expect_true(all(is.finite(b$lower) & is.finite(b$upper) & b$lower <= b$upper))
    }
  }
  b = predict(farfima(X[, 1:100], X[, 101:130], d = 0.3, order = c(1, 0), h = 3, trim = 20, band = "PID",
                       verbose = FALSE), verbose = FALSE)
  expect_equal((b$bands$PID[["80%"]]$lower + b$bands$PID[["80%"]]$upper) / 2, b$point_forecast)
})

test_that("the PID band can add a naive or theta scorecaster", {
  for (sc in c("naive", "theta")) {
    for (sym in c(TRUE, FALSE)) {
      fit = farfima(X[, 1:100], X[, 101:130], d = 0.3, order = c(1, 0), h = 3, trim = 20, band = c("NS", "PID"),
                     symmetric = sym, scorecast = sc, verbose = FALSE)
      fc = predict(fit, test = X[, 131:160], verbose = FALSE)
      b = fc$bands$PID[["80%"]]
      expect_true(all(is.finite(b$lower) & is.finite(b$upper) & b$lower <= b$upper))
      expect_true(all(is.finite(fc$evaluation$interval$MIS)))
    }
  }
})
