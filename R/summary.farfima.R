# Print the fit and the candidate table of the order selection
#' @export
summary.farfima = function(object, ...) {
  print(object)
  if (!is.null(object$order_selection)) {
    cat("\nCandidates evaluated:\n")
    print(object$order_selection$table, row.names = FALSE, digits = 4)
  }
  invisible(object)
}
