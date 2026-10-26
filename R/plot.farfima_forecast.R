#' Plot forecasts and conformal prediction bands
#'
#' * `type = "curve"`: a single plot of the forecast `h` steps beyond the data
#'   with its conformal prediction band. With `target_test_curve = k`, the
#'   forecast of the k-th test curve made `h` steps earlier, with the observed
#'   curve.
#' * `type = "horizons"`: a grid of plots, one column per band type and one row
#'   per horizon, each showing the point forecast with all levels of the band
#'   shaded. With a single band type, the horizons are side by side. With a
#'   test sample, the forecasts come from the latest test origin that reaches
#'   every horizon shown, and the observed curves are added.
#' * `type = "series"`: a row of plots, one per grid point in `at`, each showing
#'   that grid point over time with the forecasts beyond the data.
#' * `type = "recent"`: a single plot of the last `n_recent` curves in grey and
#'   the forecasts beyond the data, coloured by horizon.
#'
#' @param x A forecast from [`predict()`][predict.farfima].
#' @param type `"curve"`, `"horizons"`, `"series"` or `"recent"`.
#' @param h Forecast horizon. For `"horizons"` and `"recent"`, one or more
#'   horizons: by default 1, `h/2` and `h` for `"horizons"`, and all for
#'   `"recent"`.
#' @param band Conformal prediction band to draw. For `"horizons"`, one or
#'   more bands, by default all.
#' @param level Level of the band for `type = "curve"`.
#' @param target_test_curve Position in the test sample of the curve forecast,
#'   between `h` and the number of test curves (`type = "curve"`).
#' @param at Grid values for `type = "series"`, by default four spread over the
#'   grid.
#' @param n_recent Number of curves shown for `type = "recent"`.
#' @param argvals Grid values. Defaults to the row names of the curves if
#'   numeric.
#' @param xlab,ylab Axis labels.
#' @param main Title. For `type = "series"`, the text before the grid value in
#'   each panel title, such as `"Age"`.
#' @param legend Legend position, or `NULL` for no legend.
#' @param ... Other arguments passed to [graphics::plot()], such as `xlim`.
#' @return `x`, invisibly.
#' @seealso [`predict()`][predict.farfima]
#' @examples
#' sim = sim_farfima(n = 200, d = 0.3, p = 1, q = 0, no_grid = 51, seed = 1)
#' blocks = split_curves(sim)
#' fit = farfima(blocks$train, blocks$valid, d = 0.3, order = c(1, 0), h = 5,
#'                trim = 50, verbose = FALSE)
#' fc = predict(fit, test = blocks$test, verbose = FALSE)
#' plot(fc, h = 1)
#' plot(fc, h = 2, target_test_curve = 10)
#' plot(fc, type = "horizons")
#' plot(fc, type = "series")
#' plot(fc, type = "recent")
#' @export
plot.farfima_forecast = function(x, type = c("curve", "horizons", "series", "recent"), h = 1,
                                  band = attr(x, "settings")$band[1], level = attr(x, "settings")$level[1],
                                  target_test_curve = NULL, at = NULL, n_recent = 30, argvals = NULL,
                                  xlab = NULL, ylab = "Value", main = NULL, legend = "bottomright", ...) {
  type = match.arg(type)
  set = attr(x, "settings")
  if (length(set$band) > 0 && (!all(band %in% set$band) || (type != "horizons" && length(band) != 1))) {
    stop(sprintf("`band` must be one of the bands returned by predict(): %s.", paste(set$band, collapse = ", ")),
         call. = FALSE)
  }
  if (type == "horizons" && missing(band)) {
    band = set$band
  }
  if (type == "horizons" && missing(h)) {
    h = unique(c(1, ceiling(set$h / 2), set$h))
  }
  if (any(h < 1 | h > set$h) || (!type %in% c("horizons", "recent") && length(h) != 1)) {
    stop(sprintf("`h` must be between 1 and %d.", set$h), call. = FALSE)
  }
  if (is.null(argvals)) {
    argvals = suppressWarnings(as.numeric(rownames(x$point_forecast)))
    if (length(argvals) == 0 || any(!is.finite(argvals))) {
      argvals = seq_len(nrow(x$point_forecast))
    }
  }
  if (is.null(xlab)) {
    xlab = if (type == "series") "Time" else "Grid"
  }
  if (type == "series" && missing(legend)) {
    legend = "bottomleft"
  }

  switch(type,
         curve = plot_curve(x, h, band, level, target_test_curve, argvals, xlab, ylab, main, legend, ...),
         horizons = plot_horizons(x, h, band, argvals, xlab, ylab, main, ...),
         series = plot_series(x, band, at, argvals, xlab, ylab, main, legend, ...),
         recent = plot_recent(x, if (missing(h)) seq_len(set$h) else h, n_recent, argvals, xlab, ylab, main, legend, ...))
  invisible(x)
}

band_colours = c(SR = "#1b9e77", NS = "#d95f02", AC = "#7570b3", PID = "#e7298a", AcMCP = "#66a61e")
