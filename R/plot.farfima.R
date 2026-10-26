#' Plot the fitted or residual curves of a FARFIMA fit
#'
#' Rainbow plots of the in-sample fit of a model from [farfima()], coloured by
#' time as in [plot_curves()]. `type = "fitted"` draws the fitted curves on the
#' scale of the data, and `type = "residuals"` draws the residual curves. The
#' first `max(p, q)` curves have no fit and are left out.
#'
#' @param x A fit from [farfima()] with `method = "farma"`.
#' @param type `"fitted"` or `"residuals"`.
#' @param argvals Grid values. Defaults to the row names of the curves if
#'   numeric.
#' @param xlab,ylab,main Axis labels and title.
#' @param ... Other arguments passed to [graphics::matplot()].
#' @return `x`, invisibly.
#' @seealso [`fitted()`][fitted.farfima] and [`residuals()`][fitted.farfima]
#'   for the curves.
#' @examples
#' sim = sim_farfima(n = 200, d = 0.3, p = 1, q = 0, no_grid = 51, seed = 1)
#' fit = farfima(sim, d = 0.3, order = c(1, 0), trim = 50, verbose = FALSE)
#' plot(fit)
#' plot(fit, type = "residuals")
#' @export
plot.farfima = function(x, type = c("fitted", "residuals"), argvals = NULL, xlab = "Grid", ylab = NULL,
                        main = NULL, ...) {
  type = match.arg(type)
  curves = if (type == "fitted") fitted(x) else residuals(x)
  curves = curves[, colSums(!is.na(curves)) > 0, drop = FALSE]
  if (is.null(ylab)) {
    ylab = if (type == "fitted") "Value" else "Residual"
  }
  if (is.null(main)) {
    main = if (type == "fitted") "Fitted curves" else "Residual curves"
  }
  plot_curves(curves, argvals = argvals, xlab = xlab, ylab = ylab, main = main, ...)
  invisible(x)
}
