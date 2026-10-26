test_that("split_curves gives consecutive training, validation and test samples", {
  b = split_curves(X, train = 0.6, valid = 0.2)
  expect_named(b, c("train", "valid", "test"))
  expect_equal(c(ncol(b$train), ncol(b$valid), ncol(b$test)), c(96, 32, 32))
  expect_equal(cbind(b$train, b$valid, b$test), X, ignore_attr = TRUE)
  expect_null(split_curves(X, train = 0.8, valid = 0.2)$test)
  expect_error(split_curves(X, train = 0.8, valid = 0.4), "train \\+ valid <= 1")
})
