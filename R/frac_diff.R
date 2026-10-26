#' Fractional differencing
#'
#' Applies \eqn{(1 - B)^d} at every grid point, where \eqn{B} is the backshift
#' operator. A negative `d` integrates, so `frac_diff(frac_diff(x, d), -d)`
#' returns `x`.
#'
#' @param x A matrix of curves (grid points by time).
#' @param d Memory parameter, the order of differencing. A negative `d`
#'   fractionally integrates the curves.
#' @return The fractionally differenced curves, or the fractionally integrated
#'   curves when `d` is negative, as a matrix of the same size as `x`.
#' @seealso [estimate_memory()]
#' @examples
#' x = sim_farfima(n = 100, d = 0.3, p = 1, q = 0, seed = 1)
#' y = frac_diff(x, 0.3)
#' max(abs(frac_diff(y, -0.3) - x))
#' @export
frac_diff = function(x, d) {
  x = as.matrix(x)
  t(FCVAR::FracDiff(t(x), d = d))
}
