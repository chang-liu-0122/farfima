# Shading of each nominal level, lighter for the wider bands
level_shades = function(levels) {
  stats::setNames(seq(0.3, 0.6, length.out = length(levels)), level_label(sort(levels, decreasing = TRUE)))
}

# One forecast curve with its band, beyond the data or for a test curve
plot_curve = function(x, h, band, level, target_test_curve, argvals, xlab, ylab, main, legend, ...) {
  set = attr(x, "settings")
  origin = attr(x, "origin")
  has_band = length(set$band) > 0
  has_time = isTRUE(attr(x, "has_time"))

  if (is.null(target_test_curve)) {
    f = x$point_forecast[, h]
    b = if (has_band) x$bands[[band]][[level_label(level)]] else NULL
    ref = attr(x, "last_curve")
    ref_label = if (has_time) paste("last observed:", origin) else "last observed curve"
    ref_lty = 2
    when = if (has_time) {
      paste0(colnames(x$point_forecast)[h], " forecast from ", origin, " (h = ", h, ")")
    } else {
      paste0("h = ", h, " after the last curve")
    }
  } else {
    records = attr(x, "test_forecasts")
    if (is.null(records)) {
      stop("`target_test_curve` needs a forecast from predict() with `test`.", call. = FALSE)
    }
    n_test = attr(x, "n_test")
    k = target_test_curve
    if (length(k) != 1 || !is.numeric(k) || k < h || k > n_test) {
      stop(sprintf("With h = %d, `target_test_curve` must be between %d and %d.", h, h, n_test), call. = FALSE)
    }
    hit = Filter(function(r) length(r$target_test_curve) >= h && r$target_test_curve[h] == k, records)
    if (length(hit) == 0) {
      stop(sprintf("The %d-step forecast of test curve %d failed.", h, k), call. = FALSE)
    }
    r = hit[[1]]
    f = r$point_forecast[, h]
    b = if (has_band) r$bands[[band]][[level_label(level)]] else NULL
    ref = r$observed[, h]
    ref_label = if (has_time) paste("observed:", r$target[h]) else paste("observed test curve", k)
    ref_lty = 1
    when = if (has_time) {
      paste0(r$target[h], " forecast from ", r$origin, " (test curve ", k, ")")
    } else {
      paste0("Test curve ", k, ", forecast ", h, " step", if (h > 1) "s", " ahead")
    }
  }
  if (has_band && is.null(b)) {
    stop(sprintf("`level` must be one of %s.", paste(set$level, collapse = ", ")), call. = FALSE)
  }
  lo = if (is.null(b)) NULL else b$lower[, h]
  hi = if (is.null(b)) NULL else b$upper[, h]
  if (is.null(main)) {
    main = if (has_band) paste0(when, ", ", level_label(level), " ", band) else when
  }

  col = if (has_band) band_colours[[band]] else "steelblue4"
  graphics::plot(argvals, f, type = "n", ylim = range(c(lo, hi, f, ref), finite = TRUE),
                 xlab = xlab, ylab = ylab, main = main, ...)
  if (has_band) {
    ok = !is.na(lo) & !is.na(hi)
    graphics::polygon(c(argvals[ok], rev(argvals[ok])), c(clip_edges(lo[ok]), rev(clip_edges(hi[ok]))),
                      col = grDevices::adjustcolor(col, 0.35), border = NA)
  }
  ref_col = if (ref_lty == 1) "black" else "grey50"
  graphics::lines(argvals, ref, col = ref_col, lty = ref_lty)
  graphics::lines(argvals, f, col = col, lwd = 2)
  if (!is.null(legend)) {
    graphics::legend(legend, c(ref_label, "forecast", if (has_band) paste(level_label(level), band, "band")),
                     col = c(ref_col, col, if (has_band) grDevices::adjustcolor(col, 0.35)),
                     lty = c(ref_lty, 1, if (has_band) NA), lwd = c(1, 2, if (has_band) NA),
                     pch = c(NA, NA, if (has_band) 15), pt.cex = 2, bty = "n")
  }
}

# Clip infinite band edges to the panel
clip_edges = function(v) {
  usr = graphics::par("usr")
  pmin(pmax(v, usr[3]), usr[4])
}

# Grid of forecasts: one column per band and one row per horizon
plot_horizons = function(x, h, band, argvals, xlab, ylab, main, ...) {
  set = attr(x, "settings")
  records = Filter(function(r) length(r$target) >= max(h), attr(x, "test_forecasts"))
  if (length(records) > 0) {
    r = records[[length(records)]]
    fcst = r$point_forecast
    bands = r$bands
    observed = r$observed
    origin = r$origin
    target = r$target
  } else {
    fcst = x$point_forecast
    bands = x$bands
    observed = NULL
    origin = attr(x, "origin")
    target = colnames(fcst)
  }
  rows = if (length(set$band) > 0) band else "forecast"
  shades = level_shades(set$level)
  layout = if (length(rows) > 1) c(length(h), length(rows)) else c(1, length(h))

  if (prod(layout) > 1) {
    old = graphics::par(mfrow = layout, mar = c(4, 4, 2, 0.5), oma = c(0, 0, 2, 0))
    on.exit(graphics::par(old))
  }
  for (hh in h) {
    for (b in rows) {
      col = if (b %in% names(band_colours)) band_colours[[b]] else "steelblue4"
      edges = unlist(lapply(names(shades), function(l) c(bands[[b]][[l]]$lower[, hh], bands[[b]][[l]]$upper[, hh])))
      graphics::plot(argvals, fcst[, hh], type = "n", ylim = range(fcst[, hh], observed[, hh], edges, finite = TRUE),
                     xlab = xlab, ylab = ylab,
                     main = paste0(if (b %in% names(band_colours)) paste0(b, ": "), target[hh], " (h = ", hh, ")"), ...)
      for (l in names(shades)) {
        lo = bands[[b]][[l]]$lower[, hh]
        up = bands[[b]][[l]]$upper[, hh]
        if (!is.null(lo)) {
          graphics::polygon(c(argvals, rev(argvals)), c(clip_edges(lo), rev(clip_edges(up))),
                            col = grDevices::adjustcolor(col, shades[[l]]), border = NA)
          graphics::lines(argvals, lo, col = col, lwd = 0.8)
          graphics::lines(argvals, up, col = col, lwd = 0.8)
        }
      }
      graphics::lines(argvals, fcst[, hh], col = "grey20", lwd = 1.5, lty = 2)
      if (!is.null(observed)) {
        graphics::lines(argvals, observed[, hh], lwd = 1.5)
      }
    }
  }
  if (is.null(main)) {
    main = paste0("Forecasts from ", origin, if (length(set$level) > 0) paste0(", ", paste(level_label(set$level), collapse = " and "), " bands"),
                  if (!is.null(observed)) "; black: observed", "; dashed: point forecast")
  }
  if (prod(layout) > 1) {
    graphics::mtext(main, outer = TRUE, line = 0.5, cex = 0.9)
  } else {
    graphics::mtext(main, side = 3, line = 2.5, cex = 0.9)
  }
}

# Curves over time at a few grid points
plot_series = function(x, band, at, argvals, xlab, ylab, main, legend, ...) {
  set = attr(x, "settings")
  X = attr(x, "curves")
  times = suppressWarnings(as.numeric(colnames(X)))
  if (length(times) == 0 || any(!is.finite(times))) {
    times = seq_len(ncol(X))
  }
  future = suppressWarnings(as.numeric(colnames(x$point_forecast)))
  if (any(!is.finite(future))) {
    future = times[length(times)] + seq_len(ncol(x$point_forecast))
  }
  if (is.null(at)) {
    idx = unique(round(seq(1, nrow(X), length.out = 4)))
  } else {
    idx = match(at, argvals)
    if (any(is.na(idx))) {
      stop("`at` must be values of the grid.", call. = FALSE)
    }
  }
  has_band = length(set$band) > 0
  col = if (has_band) band_colours[[band]] else "steelblue4"
  shades = level_shades(set$level)

  if (length(idx) > 1) {
    old = graphics::par(mfrow = c(ceiling(length(idx) / 4), min(length(idx), 4)), mar = c(4, 4, 2.5, 0.5))
    on.exit(graphics::par(old))
  }
  args = list(...)
  xlim = if (is.null(args$xlim)) range(times, future) else args$xlim
  shown = times >= xlim[1] & times <= xlim[2]
  args$xlim = NULL
  for (k in idx) {
    edges = if (has_band) unlist(lapply(x$bands[[band]], function(b) c(b$lower[k, ], b$upper[k, ]))) else NULL
    ylim = if (is.null(args$ylim)) range(X[k, shown], x$point_forecast[k, ], edges, finite = TRUE) else args$ylim
    do.call(graphics::plot, c(list(times, X[k, ], type = "l", xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab,
                                   main = paste(if (is.null(main)) "Grid point" else main, argvals[k])),
                              args[names(args) != "ylim"]))
    t_last = times[length(times)]
    y_last = X[k, ncol(X)]
    if (has_band) {
      for (l in names(shades)) {
        b = x$bands[[band]][[l]]
        graphics::polygon(c(t_last, future, rev(future)), c(y_last, clip_edges(b$lower[k, ]), rev(clip_edges(b$upper[k, ]))),
                          col = grDevices::adjustcolor(col, shades[[l]]), border = NA)
      }
    }
    graphics::lines(c(t_last, future), c(y_last, x$point_forecast[k, ]), col = col, lwd = 2)
  }
  if (!is.null(legend) && has_band) {
    graphics::legend(legend, c("observed", "forecast", paste(names(shades), band, "band")),
                     col = c("black", col, sapply(shades, function(a) grDevices::adjustcolor(col, a))),
                     lty = c(1, 1, rep(NA, length(shades))),
                     lwd = c(1, 2, rep(NA, length(shades))), pch = c(NA, NA, rep(15, length(shades))), pt.cex = 2,
                     bty = "n", cex = 0.8)
  }
}

# Recent curves and the forecasts beyond the data
plot_recent = function(x, h, n_recent, argvals, xlab, ylab, main, legend, ...) {
  h = sort(unique(h))
  X = attr(x, "curves")
  recent = X[, max(1, ncol(X) - n_recent + 1):ncol(X), drop = FALSE]
  fcst = x$point_forecast[, h, drop = FALSE]
  cols = grDevices::hcl.colors(ncol(x$point_forecast), "Plasma")[h]
  if (is.null(main)) {
    which_h = if (length(h) == 1) paste("h =", h) else if (all(diff(h) == 1)) paste0("h = ", h[1], " to ", max(h)) else paste("h =", paste(h, collapse = ", "))
    main = paste0("Last ", ncol(recent), " observed curves and the forecast", if (length(h) > 1) "s", " for ", which_h)
  }
  graphics::plot(argvals, fcst[, 1], type = "n", ylim = range(recent, fcst, finite = TRUE),
                 xlab = xlab, ylab = ylab, main = main, ...)
  graphics::matlines(argvals, recent, lty = 1, col = "grey75")
  graphics::matlines(argvals, fcst, lty = 1, lwd = 1.5, col = cols)
  if (!is.null(legend)) {
    ends = if (length(h) <= 4) seq_along(h) else c(1, length(h))
    graphics::legend(legend, c("observed", paste0("h = ", h[ends], " (", colnames(fcst)[ends], ")")),
                     col = c("grey60", cols[ends]), lty = 1, lwd = c(1, rep(2, length(ends))), bty = "n")
  }
}
