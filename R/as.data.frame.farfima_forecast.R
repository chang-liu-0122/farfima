#' Forecasts and conformal prediction bands as a data frame
#'
#' One row per grid point, horizon, band and level, for saving the results or
#' for your own plots and summaries.
#'
#' @param x A forecast from [`predict()`][predict.farfima].
#' @param row.names,optional Not used; they come from [as.data.frame()].
#' @param test `FALSE` for the forecasts after the last validation or test
#'   curve; `TRUE` for the forecasts of the test curves, with the observed
#'   values.
#' @param ... Not used.
#' @return A data frame with one row per grid point, horizon, band and level:
#'   * `band`, `level`, `h`: the band type, its nominal level and the horizon;
#'   * `target`: the time being forecast;
#'   * `grid`: the grid point;
#'   * `point_forecast`, `lower`, `upper`: the point forecast and the edges of
#'     the conformal prediction band at that grid point.
#'
#'   With `test = TRUE`, the rows cover the forecasts of the test curves, with
#'   three more columns:
#'   * `origin`: the time the forecast was made from;
#'   * `target_test_curve`: the position of the target in the test sample;
#'   * `observed`: the observed value at that grid point.
#' @seealso [`predict()`][predict.farfima]
#' @examples
#' sim = sim_farfima(n = 200, d = 0.3, p = 1, q = 0, no_grid = 51, seed = 1)
#' blocks = split_curves(sim)
#' fit = farfima(blocks$train, blocks$valid, d = 0.3, order = c(1, 0), h = 3,
#'                trim = 50, verbose = FALSE)
#' fc = predict(fit, test = blocks$test, verbose = FALSE)
#' head(as.data.frame(fc))
#' d = as.data.frame(fc, test = TRUE)
#' d = subset(d, band == "AC" & level == 0.8)
#' tapply(d$observed >= d$lower & d$observed <= d$upper, d$h, mean)
#' @export
as.data.frame.farfima_forecast = function(x, row.names = NULL, optional = FALSE, test = FALSE, ...) {
  grid = if (is.null(rownames(x$point_forecast))) seq_len(nrow(x$point_forecast)) else rownames(x$point_forecast)

  # One long table for a set of forecasts and their bands
  long = function(pf, bands, target) {
    base = data.frame(h = rep(seq_along(target), each = length(grid)), target = rep(target, each = length(grid)),
                      grid = rep(grid, length(target)), point_forecast = as.vector(pf))
    if (length(bands) == 0) {
      return(base)
    }
    do.call(rbind, lapply(names(bands), function(m) {
      do.call(rbind, lapply(names(bands[[m]]), function(lev) {
        b = bands[[m]][[lev]]
        data.frame(band = m, level = as.numeric(sub("%", "", lev)) / 100, base,
                   lower = as.vector(b$lower), upper = as.vector(b$upper))
      }))
    }))
  }

  if (!test) {
    return(long(x$point_forecast, x$bands, colnames(x$point_forecast)))
  }
  records = attr(x, "test_forecasts")
  if (is.null(records)) {
    stop("This forecast has no test evaluation.", call. = FALSE)
  }
  do.call(rbind, lapply(records, function(r) {
    out = data.frame(origin = r$origin, long(r$point_forecast, r$bands, r$target))
    out$target_test_curve = r$target_test_curve[out$h]
    out$observed = rep(as.vector(r$observed), length.out = nrow(out))
    out
  }))
}
