test_that("sim_farfima gives curves of the requested size and is reproducible with a seed", {
  x = sim_farfima(n = 50, d = 0.2, p = 1, q = 1, no_grid = 11, seed = 2)
  expect_equal(dim(x), c(11, 50))
  expect_equal(dim(attr(x, "innovations")), c(11, 50))
  expect_identical(x, sim_farfima(n = 50, d = 0.2, p = 1, q = 1, no_grid = 11, seed = 2))
  expect_error(sim_farfima(n = 50, d = 0.2, p = 0, q = 0), "cannot both be zero")
})
