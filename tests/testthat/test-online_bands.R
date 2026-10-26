test_that("online_bands uses only the curves observed by each origin", {
  set.seed(3)
  F = lapply(1:2, function(hh) matrix(rnorm(21 * 12), 21))
  Y = lapply(1:2, function(hh) F[[hh]] + matrix(rnorm(21 * 12), 21))
  O = list(1:12, 1:12)
  ref = online_reference(lapply(1:2, function(hh) matrix(rnorm(21 * 30), 21)), list(-29:0, -29:0))
  for (m in c("AC", "PID", "AcMCP")) {
    a = online_bands(m, FALSE, F, Y, O, ref, 0.2, 0.1)
    Y2 = Y
    Y2[[2]][, 6] = Y2[[2]][, 6] + 10
    b = online_bands(m, FALSE, F, Y2, O, ref, 0.2, 0.1)
    expect_equal(a[[2]]$lower[, 1:7], b[[2]]$lower[, 1:7])
    expect_false(isTRUE(all.equal(a[[2]]$lower[, 8], b[[2]]$lower[, 8])))
    longer = online_bands(m, FALSE, lapply(F, cbind, 0), lapply(Y, cbind, 1), list(1:13, 1:13), ref, 0.2, 0.1)
    expect_equal(longer[[2]]$lower[, 1:12], a[[2]]$lower)
  }
})

test_that("the symmetric PID band without learning or integral term is the NS band", {
  set.seed(3)
  F = lapply(1:2, function(hh) matrix(rnorm(21 * 12), 21))
  Y = lapply(1:2, function(hh) F[[hh]] + matrix(rnorm(21 * 12), 21))
  O = list(1:12, 1:12)
  ref = online_reference(lapply(1:2, function(hh) matrix(rnorm(21 * 30), 21)), list(-29:0, -29:0))
  p = online_bands("PID", TRUE, F, Y, O, ref, 0.2, 0, integrate = FALSE)
  q = capped_quantile(sort(abs(ref$h[[1]]$z)), 0.2)
  expect_equal(p[[1]]$upper - F[[1]], q * matrix(ref$h[[1]]$gamma, 21, 12))
  small = online_bands("PID", TRUE, lapply(F, function(f) f[1:3, ]), lapply(Y, function(y) y[1:3, ]), O,
                       online_reference(lapply(1:2, function(hh) matrix(rnorm(3 * 30), 3)), list(-29:0, -29:0)), 0.01, 0.1)
  expect_true(all(is.finite(small[[2]]$upper)))
})

test_that("pid_integrator saturates to an infinite band", {
  expect_equal(pid_integrator(10, 20, 1), Inf)
  expect_equal(pid_integrator(-10, 20, 1), -Inf)
  expect_equal(pid_integrator(1, 20, 2), 2 * tan(log(20) / (2 / pi * (1 - 1 / log(20)) * 20)))
})

test_that("the AcMCP centre corrects a persistent bias", {
  set.seed(4)
  u = seq(0, 1, length.out = 21)
  make = function(n) {
    F = lapply(1:2, function(hh) matrix(rnorm(21 * n), 21))
    Y = lapply(1:2, function(hh) F[[hh]] + hh * sin(pi * u) + matrix(rnorm(21 * n, sd = 0.3), 21))
    list(F = F, Y = Y, E = lapply(1:2, function(hh) Y[[hh]] - F[[hh]]))
  }
  cal = make(40)
  run = make(60)
  ref = online_reference(cal$E, list(1:40, 1:40))
  out = online_bands("AcMCP", FALSE, run$F, run$Y, list(41:100, 41:100), ref, 0.2, 0.05)
  ecp = sapply(out, function(r) r$tally[1] / r$tally[3])
  expect_true(all(abs(ecp - 0.8) < 0.05))
  expect_equal(rowMeans(out[[2]]$centre * ref$h[[2]]$gamma), 2 * sin(pi * u), tolerance = 0.1)
  sym = online_bands("PID", TRUE, run$F, run$Y, list(41:100, 41:100), ref, 0.2, 0.05)
  expect_lt(out[[2]]$tally[2], sym[[2]]$tally[2])
})

test_that("score_forecast gives the naive and theta scorecasts", {
  Z = matrix(c(1, 2, 3, 4, 5, 2, 2, 2, 2, 2), 2, byrow = TRUE)
  expect_equal(score_forecast(Z, 2, "naive"), c(5, 2))
  expect_equal(score_forecast(Z[, 1:2], 2, "theta"), c(2, 2))
  f = score_forecast(Z, 1, "theta")
  expect_equal(f[2], 2)
  expect_true(f[1] > 5 && f[1] < 6)
})

test_that("ma_intercept is the intercept of an MA(h - 1) fit", {
  x = c(0.5, -0.2, 0.7, 0.4, 0.8, -0.3, 0.6, -0.4)
  expect_equal(ma_intercept(x, 1), mean(x))
  expect_equal(ma_intercept(x, 2), unname(stats::arima(x, order = c(0, 0, 1), method = "CSS-ML")$coef[["intercept"]]))
  expect_true(is.na(ma_intercept(x[1:3], 3)))
  set.seed(5)
  y = 1 + stats::arima.sim(list(ma = c(0.5, 0.3)), n = 300)
  expect_equal(ma_intercept(y, 3), 1, tolerance = 0.2)
})
