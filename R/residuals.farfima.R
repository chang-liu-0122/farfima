# Residual curves of the fit
#' @rdname fitted.farfima
#' @export
residuals.farfima = function(object, ...) {
  if (attr(object, "settings")$method == "far") {
    stop("Residuals are only available for the FARMA method.", call. = FALSE)
  }
  object$residuals
}
