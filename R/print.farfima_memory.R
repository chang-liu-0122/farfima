# Print a memory estimate
#' @export
print.farfima_memory = function(x, ...) {
  cat("Memory estimate d = ", format(x$d, digits = 4), "\n", sep = "")
  cat("  ", memory_label(x), "\n", sep = "")
  if (isTRUE(is.finite(x$bias))) {
    cat("  uncorrected d = ", format(x$d_initial, digits = 4),
        ", bootstrap bias = ", format(x$bias, digits = 3),
        ", bootstrap sd = ", format(attr(x, "sd"), digits = 3), "\n", sep = "")
  }
  if (!isTRUE(attr(x, "diagnostics")$admissible)) {
    cat("  inadmissible: the bias-corrected estimate is below zero\n")
  }
  cat(memory_notes(x, "  "), sep = "")
  invisible(x)
}
