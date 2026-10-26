# Print the selected orders, with a note on candidates that failed during CV
#' @export
print.farfima_order = function(x, ...) {
  cat("Selected (p, q) = (", x$p, ", ", x$q, ")\n", sep = "")
  cat("  ", order_label(x), "\n", sep = "")
  if (attr(x, "method") == "cv") {
    failed = sum(!is.finite(x$table$RMSPE))
    if (failed > 0) {
      cat("  note: ", failed, " of ", nrow(x$table), " candidates failed at one or more validation origins and were dropped\n",
          sep = "")
    }
  }
  invisible(x)
}
