# Print the forecast settings and the test evaluation averaged over horizons
#' @export
print.farfima_forecast = function(x, ...) {
  m = attr(x, "model")
  set = attr(x, "settings")
  cat("Forecasts from FARFIMA(p = ", m$order[["p"]], ", d = ", format(m$d, digits = 4), ", q = ", m$order[["q"]],
      ")\n", sep = "")
  cat("  Horizons:  h = 1..", set$h, " after curve ", attr(x, "origin"), " (", set$window, " window)\n", sep = "")
  if (length(set$band) > 0) {
    cat("  Bands:     ", paste0(set$band, " (", band_names[set$band], ")", collapse = ", "), "\n", sep = "")
    cat("  Levels:    ", paste(level_label(set$level), collapse = ", "), "\n", sep = "")
  } else {
    cat("  Bands:     none (the fit has no validation sample)\n")
  }

  if (!is.null(x$evaluation)) {
    ev = x$evaluation
    cat("\nTest evaluation over ", attr(x, "n_test"), " test curves, averaged over horizons:\n", sep = "")
    cat("  RMSPE = ", format(mean(ev$point$RMSPE), digits = 4),
        ", MAPE = ", format(mean(ev$point$MAPE), digits = 4), "\n", sep = "")
    if (!is.null(ev$interval)) {
      tab = stats::aggregate(ev$interval[, c("ECP", "CPD", "MIS")], by = ev$interval[, c("band", "level")], FUN = mean)
      tab = tab[order(tab$level, match(tab$band, c("SR", "NS", "AC", "PID", "AcMCP"))), ]
      tab$level = level_label(tab$level)
      print(tab, row.names = FALSE, digits = 4)
      w = tapply(ev$wis$WIS, ev$wis$band, mean)[set$band]
      cat("  WIS: ", paste0(names(w), " = ", format(w, digits = 4), collapse = ", "), "\n", sep = "")
    }
  }
  invisible(x)
}
