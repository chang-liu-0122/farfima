# Print the forecast and the point accuracy by horizon
#' @export
summary.farfima_forecast = function(object, ...) {
  print(object)
  if (!is.null(object$evaluation)) {
    cat("\nPoint accuracy by horizon:\n")
    print(object$evaluation$point, row.names = FALSE, digits = 4)
  }
  invisible(object)
}
