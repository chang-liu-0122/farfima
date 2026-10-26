# Print a fitted FARFIMA model
#' @export
print.farfima = function(x, ...) {
  set = attr(x, "settings")
  dat = attr(x, "data")
  cat("FARFIMA(p = ", x$order[["p"]], ", d = ", format(x$d, digits = 4), ", q = ", x$order[["q"]], ")\n", sep = "")
  cat("  Data:      ", nrow(dat$train), " grid points x ", ncol(dat$train), " training curves",
      if (!is.null(dat$valid)) paste0(", ", ncol(dat$valid), " validation curves"), "\n", sep = "")

  cat("  Memory:    ", if (is.null(x$memory)) "supplied" else memory_label(x$memory), "\n", sep = "")
  if (isTRUE(is.finite(x$memory$bias))) {
    cat("             uncorrected d = ", format(x$memory$d_initial, digits = 4),
        ", bootstrap bias = ", format(x$memory$bias, digits = 3), "\n", sep = "")
  }
  if (!is.null(x$memory)) {
    cat(memory_notes(x$memory, "             "), sep = "")
  }

  if (set$method == "far") {
    orders = "fixed at (1, 0) for the FAR(1) benchmark"
    operators = paste0("score-based FAR(1) with ", attr(x, "estimator")$kn, " components (",
                       100 * set$far_pve, "% of the variation)")
  } else {
    orders = if (is.null(x$order_selection)) "supplied" else order_label(x$order_selection)
    operators = paste0("functional Hannan-Rissanen estimator with ", set$Npc, " principal component",
                       if (set$Npc > 1) "s")
  }
  cat("  Orders:    ", orders, "\n", sep = "")
  cat("  Operators: ", operators, ", first ", set$trim, " fractionally differenced curves trimmed\n", sep = "")

  reasons = attr(x, "stability")$reasons
  if (length(reasons) > 0) {
    cat("  Stability: unstable (", paste(reasons, collapse = ", "), ")\n", sep = "")
  }

  if (is.null(x$calibration)) {
    cat("  Bands:     none (no validation sample), point forecasts only\n")
  } else {
    asym = intersect(set$band, c("AC", "PID"))
    notes = c(if (isFALSE(set$symmetric) && length(asym) > 0) paste(paste(asym, collapse = " and "), "asymmetric"),
              if ("PID" %in% set$band && !is.null(set$scorecast) && set$scorecast != "none")
                paste0("PID with ", set$scorecast, " scorecaster"))
    cat("  Bands:     ", paste(set$band, collapse = ", "),
        if (length(notes) > 0) paste0(" (", paste(notes, collapse = "; "), ")"),
        " at ", paste(level_label(set$level), collapse = ", "),
        " for h = 1..", set$h, ", calibrated on ", ncol(dat$valid), " validation curves (", set$window,
        " window)\n", sep = "")
    n_failed = length(attr(x, "failures"))
    if (n_failed > 0) {
      cat("             ", n_failed, " of ", ncol(dat$valid), " validation forecasts failed and were skipped\n", sep = "")
    }
  }
  invisible(x)
}
