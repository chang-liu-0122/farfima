# Print the orders and the operator norm of each operator
#' @export
print.farfima_operators = function(x, ...) {
  ops = operator_list(x)
  Npc = attr(x, "Npc")
  cat("FARMA(", ops$p, ", ", ops$q, ") operators\n", sep = "")
  cat("  Data:      ", nrow(x$residuals), " grid points x ", ncol(x$residuals), " curves\n", sep = "")
  cat("  Operators: functional Hannan-Rissanen estimator with ", Npc, " principal component",
      if (Npc > 1) "s", "\n", sep = "")
  norms = vapply(ops$matrices, function(A) format(norm(A, type = "2"), digits = 3), character(1))
  cat("  Norms:     ", paste(sub("operator, lag", "lag", ops$titles), norms, sep = " = ", collapse = ", "), "\n",
      sep = "")
  reasons = attr(x, "stability")$reasons
  if (length(reasons) > 0) {
    cat("  Stability: unstable (", paste(reasons, collapse = ", "), ")\n", sep = "")
  }
  invisible(x)
}
