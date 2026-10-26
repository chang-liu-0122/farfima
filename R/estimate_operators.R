#' Estimate the FARMA operators
#'
#' Estimates the operators of a FARMA(p, q) model by the functional
#' Hannan-Rissanen method of Kuenzer (2024). The curves are assumed to be
#' short-memory, so long-memory curves need to be fractionally differenced
#' first, for example with [frac_diff()]. [farfima()] does this on the centred
#' training curves and drops the first `trim` before estimating the operators.
#'
#' Each operator is stored as a grid by grid matrix, the kernel times the grid
#' spacing, so that applying the operator to a curve is a matrix product.
#' `plot()` draws the operators as heatmaps, red for positive and blue for
#' negative values.
#'
#' The estimation code is adapted from code by Thomas Kuenzer.
#'
#' @param x A matrix of short-memory curves (grid points by time) or an `fts`
#'   object. For `plot()`, a `farfima_operators` object.
#' @param order FARMA orders `c(p, q)`, or a result of [select_order()].
#' @param Npc Number of principal components in the dimension reduction of the
#'   Hannan-Rissanen method.
#' @param argvals Grid values. Defaults to the row names of the curves if
#'   numeric.
#' @param xlab,ylab,main Axis labels and title.
#' @param ... Other arguments passed to [graphics::image()].
#' @return An object of class `farfima_operators` with `order`, the AR and MA
#'   operators `phi` and `theta` (grid by grid by lag, `NULL` when the order is
#'   zero), and the one-step-ahead in-sample `fitted` curves and `residuals`.
#'   Printing it shows the norm of each operator.
#' @seealso [select_order()], [farfima()]
#' @references
#' Kuenzer, T. (2024). Estimation of functional ARMA models. *Bernoulli*.
#' @examples
#' x = sim_farfima(n = 300, d = 0, p = 1, q = 1, seed = 1)
#' op = estimate_operators(x, order = c(1, 1))
#' op
#' plot(op)
#' @export
estimate_operators = function(x, order, Npc = 1) {
  x = as_curves(x, "x")
  if (inherits(order, "farfima_order")) {
    order = c(order$p, order$q)
  }
  if (!is.numeric(order) || length(order) != 2 || any(order < 0) || all(order == 0)) {
    stop("`order` must be c(p, q) with p + q > 0, or a `farfima_order` object.", call. = FALSE)
  }
  operators = operator_fit(x, as.integer(order[1]), as.integer(order[2]), Npc)
  if (length(attr(operators, "stability")$reasons) > 0) {
    warning(sprintf("The fitted FARMA(%d, %d) is unstable; try smaller p and q.", order[1], order[2]),
            call. = FALSE)
  }
  operators
}
