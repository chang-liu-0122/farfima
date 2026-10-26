# Memory parameter and operators of the fit
#' @export
coef.farfima = function(object, ...) {
  list(d = object$d, phi = object$phi, theta = object$theta)
}
