# Heatmaps of the operators
#' @rdname estimate_operators
#' @export
plot.farfima_operators = function(x, argvals = NULL, xlab = "Grid point v", ylab = "Grid point u", main = NULL, ...) {
  ops = operator_list(x)
  if (is.null(argvals)) {
    argvals = suppressWarnings(as.numeric(rownames(x$residuals)))
    if (length(argvals) != nrow(x$residuals) || any(!is.finite(argvals))) {
      argvals = seq_len(nrow(x$residuals))
    }
  }
  limit = max(abs(unlist(ops$matrices)))
  breaks = seq(-limit, limit, length.out = 102)
  cols = grDevices::hcl.colors(101, "Blue-Red 3")
  old = graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old))

  # Up to three panels in a row, otherwise AR and MA in separate rows
  n = ops$p + ops$q
  if (n <= 3 || ops$p == 0 || ops$q == 0) {
    panel = matrix(seq_len(n), 1)
  } else {
    panel = matrix(0, 2, max(ops$p, ops$q))
    panel[1, seq_len(ops$p)] = seq_len(ops$p)
    panel[2, seq_len(ops$q)] = ops$p + seq_len(ops$q)
  }
  graphics::layout(cbind(panel, n + 1), widths = c(rep(1, ncol(panel)), graphics::lcm(2.4)))
  graphics::par(mar = c(1.2, 1.2, 2, 0.6), mgp = c(2, 0.6, 0),
                oma = c(3.5, 3.5, if (is.null(main)) 0.5 else 2.5, 0.5))
  for (k in seq_len(n)) {
    at = which(panel == k, arr.ind = TRUE)
    below = at[1] < nrow(panel) && panel[at[1] + 1, at[2]] > 0
    graphics::image(argvals, argvals, t(ops$matrices[[k]]), col = cols, breaks = breaks, axes = FALSE,
                    xlab = "", ylab = "", main = ops$titles[k], ...)
    graphics::box()
    graphics::axis(1, labels = !below)
    graphics::axis(2, labels = at[2] == 1, las = 1)
  }

  graphics::par(mar = c(1.2, 0.3, 2, 2.8))
  mids = (breaks[-1] + breaks[-length(breaks)]) / 2
  graphics::image(1, mids, matrix(mids, 1), col = cols, breaks = breaks, axes = FALSE, xlab = "", ylab = "",
                  main = "Value", cex.main = 0.9, font.main = 1)
  graphics::axis(4, las = 1)
  graphics::box()

  graphics::mtext(xlab, side = 1, outer = TRUE, line = 2)
  graphics::mtext(ylab, side = 2, outer = TRUE, line = 2)
  if (!is.null(main)) {
    graphics::mtext(main, side = 3, outer = TRUE, line = 0.8, font = 2)
  }
  invisible(x)
}
