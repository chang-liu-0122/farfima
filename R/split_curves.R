#' Split curves into training, validation and test samples
#'
#' Splits the curves in time order.
#'
#' @param x A matrix of curves (grid points by time) or an `fts` object.
#' @param train,valid Proportions of the curves in the training and validation
#'   samples, rounded down. The remaining curves form the test sample.
#' @return A list with `train`, `valid` and `test`. `test` is `NULL` when
#'   `train + valid = 1`.
#' @seealso [farfima()]
#' @examples
#' x = sim_farfima(n = 100, d = 0.2, p = 1, q = 0, no_grid = 21, seed = 1)
#' blocks = split_curves(x, train = 0.6, valid = 0.2)
#' sapply(blocks, ncol)
#' @export
split_curves = function(x, train = 0.6, valid = 0.2) {
  x = as_curves(x, "x")
  if (train <= 0 || valid <= 0 || train + valid > 1) {
    stop("`train` and `valid` must be positive proportions with train + valid <= 1.", call. = FALSE)
  }

  train_end = floor(ncol(x) * train)
  valid_end = train_end + floor(ncol(x) * valid)

  list(
    train = x[, seq_len(train_end), drop = FALSE],
    valid = x[, (train_end + 1):valid_end, drop = FALSE],
    test = if (valid_end < ncol(x)) x[, (valid_end + 1):ncol(x), drop = FALSE] else NULL
  )
}
