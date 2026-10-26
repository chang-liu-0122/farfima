test_that("as.data.frame lays out the forecasts and bands as one long table", {
  fit = farfima(X[, 1:100], X[, 101:130], d = 0.3, order = c(1, 0), h = 3, trim = 20, verbose = FALSE)
  fc = predict(fit, test = X[, 131:160], verbose = FALSE)
  expect_equal(nrow(as.data.frame(fc, test = TRUE)), (28 * 3 + 2 + 1) * 21 * 3 * 2)
  for (sym in c(TRUE, FALSE)) {
    fit = farfima(X[, 1:100], X[, 101:130], d = 0.3, order = c(1, 0), h = 3, trim = 20,
                   band = c("PID", "AcMCP"), symmetric = sym, verbose = FALSE)
    d = as.data.frame(predict(fit, test = X[, 131:160], verbose = FALSE), test = TRUE)
    expect_false(anyNA(d$lower[d$band %in% c("PID", "AcMCP")]))
  }
})
