test_that("order selection stays inside the candidate set", {
  Y = frac_diff(X - rowMeans(X), d = 0.3)[, -(1:20)]
  for (m in c("ic_full", "ic_auto")) {
    o = select_order(Y, method = m, p = 1:2, q = 0:1)
    expect_true(o$p %in% 1:2 && o$q %in% 0:1)
    expect_true(is.integer(o$p) && is.integer(o$q))
  }
  Z = frac_diff(X[, 1:130] - rowMeans(X[, 1:100]), d = 0.3)[, -(1:20)]
  cv = select_order(Z[, 1:80], method = "cv", p = 0:1, q = 0:1, valid = Z[, 81:110], h = 3)
  expect_equal(nrow(cv$table), 3)
  expect_named(cv, c("p", "q", "table"))
  grDevices::pdf(NULL)
  expect_silent(plot(cv))
  expect_silent(plot(select_order(Y, method = "ic_auto", p = 1:2, q = 0:1)))
  expect_silent(plot(select_order(Y, method = "bj", p = 0:1, q = 0:1, seed = 1)))
  grDevices::dev.off()
  expect_named(select_order(Y, method = "ic_auto", p = 1:2, q = 0:1), c("p", "q", "table", "search_path"))
})

test_that("farfima selects the orders of the fractionally differenced curves", {
  fit = farfima(X[, 1:100], d = 0.3, p = 0:1, q = 0:1, trim = 20, verbose = FALSE)
  Y = frac_diff(X[, 1:100] - rowMeans(X[, 1:100]), d = 0.3)[, -(1:20)]
  o = select_order(Y, p = 0:1, q = 0:1)
  expect_equal(fit$order, c(p = o$p, q = o$q))
  expect_equal(fit$order_selection$table, o$table)
})

test_that("cross-validation reports candidates that failed", {
  real = forecast_original_scale
  local_mocked_bindings(forecast_original_scale = function(..., q_est) {
    if (q_est == 1) stop(unstable(1.3))
    real(..., q_est = q_est)
  })
  Z = frac_diff(X[, 1:130] - rowMeans(X[, 1:100]), d = 0.3)[, -(1:20)]
  cv = select_order(Z[, 1:80], method = "cv", p = 0:1, q = 0:1, valid = Z[, 81:110], h = 3)
  expect_equal(c(cv$p, cv$q), c(1L, 0L))
  expect_output(print(cv), "2 of 3 candidates failed at one or more validation origins")
})
