#' Rainbow plots of curves
#'
#' Draws the curves coloured by time, from red (earliest) to violet (latest).
#' `type = "rainbow"` draws all curves in one panel with a colour bar, and
#' `type = "periods"` draws one panel per period with the other curves in grey.
#'
#' @param x A matrix of curves (grid points by time) or an `fts` object.
#' @param type `"rainbow"` or `"periods"`.
#' @param periods Number of periods of equal length, or the first time point
#'   of each period, such as `c(1901, 1921, 1941, 1961, 1981)`.
#' @param argvals Grid values. Defaults to the row names of `x` if numeric.
#' @param xlab,ylab,main Axis labels and title.
#' @param ... Other arguments passed to [graphics::matplot()].
#' @return `x`, invisibly.
#' @seealso [split_curves()]
#' @examples
#' x = sim_farfima(n = 100, d = 0.3, p = 1, q = 0, no_grid = 21, seed = 1)
#' plot_curves(x)
#' plot_curves(x, type = "periods", periods = 4)
#' @export
plot_curves = function(x, type = c("rainbow", "periods"), periods = 5, argvals = NULL,
                       xlab = "Grid", ylab = "Value", main = NULL, ...) {
  type = match.arg(type)
  x = as_curves(x, "x")
  n = ncol(x)
  times = suppressWarnings(as.numeric(colnames(x)))
  if (length(times) != n || any(!is.finite(times))) {
    times = seq_len(n)
  }
  if (is.null(argvals)) {
    argvals = suppressWarnings(as.numeric(rownames(x)))
    if (length(argvals) != nrow(x) || any(!is.finite(argvals))) {
      argvals = seq_len(nrow(x))
    }
  }
  cols = grDevices::rainbow(n, end = 0.8)
  old = graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old))

  # All curves, with a colour bar of the time points
  if (type == "rainbow") {
    graphics::layout(matrix(1:2, 1), widths = c(6, 1))
    graphics::par(mar = c(4, 4, 3, 0.5))
    graphics::matplot(argvals, x, type = "l", lty = 1, col = cols, xlab = xlab, ylab = ylab,
                      main = if (is.null(main)) paste0("Curves from ", times[1], " to ", times[n]) else main, ...)
    graphics::par(mar = c(4, 0.5, 3, 3.5))
    graphics::image(1, times, matrix(times, 1), col = cols, axes = FALSE, xlab = "", ylab = "")
    ticks = seq(times[1], times[n], length.out = min(n, 6))
    if (all(times == round(times))) {
      ticks = round(ticks)
    }
    graphics::axis(4, at = ticks, labels = ticks, las = 1)
    return(invisible(x))
  }

  # One panel per period
  if (length(periods) == 1) {
    group = ceiling(seq_len(n) * periods / n)
  } else {
    group = pmax(findInterval(times, sort(periods)), 1)
  }
  panels = sort(unique(group))
  n_col = min(length(panels), 3)
  graphics::par(mfrow = c(ceiling(length(panels) / n_col), n_col), mar = c(1.2, 1.2, 2, 0.6), mgp = c(2, 0.6, 0),
                oma = c(3.5, 3.5, if (is.null(main)) 0.5 else 2.5, 0.5))
  for (k in seq_along(panels)) {
    sel = group == panels[k]
    graphics::matplot(argvals, x, type = "n", xlab = "", ylab = "", axes = FALSE,
                      main = paste0(min(times[sel]), "-", max(times[sel])), ...)
    graphics::box()
    graphics::axis(1, labels = k + n_col > length(panels))
    graphics::axis(2, labels = (k - 1) %% n_col == 0, las = 1)
    graphics::matlines(argvals, x[, !sel, drop = FALSE], lty = 1, col = "grey88")
    graphics::matlines(argvals, x[, sel, drop = FALSE], lty = 1, col = cols[sel])
  }
  graphics::mtext(xlab, side = 1, outer = TRUE, line = 2)
  graphics::mtext(ylab, side = 2, outer = TRUE, line = 2)
  if (!is.null(main)) {
    graphics::mtext(main, side = 3, outer = TRUE, line = 0.8, font = 2)
  }
  invisible(x)
}
