#' Fitted curves and residuals of a FARFIMA fit
#'
#' The in-sample fit of a model from [farfima()]. Each fitted curve is the
#' one-step-ahead forecast of a training curve from the curves before it, and
#' each residual curve is the training curve minus its fitted curve. They cover
#' the training curves left after trimming; the first `max(p, q)` columns are
#' `NA`, because the model needs that many earlier curves.
#'
#' `fitted(fit)` is on the scale of the data, so `fitted(fit) + residuals(fit)`
#' gives the training curves. `fitted(fit, scale = "differenced")` is the FARMA
#' fit of the fractionally differenced curves. The residuals are the same on
#' both scales.
#'
#' @param object A fit from [farfima()] with `method = "farma"`.
#' @param scale `"original"` or `"differenced"`.
#' @param ... Not used.
#' @return `fitted()` returns the fitted curves and `residuals()` the residual
#'   curves, each as a matrix with grid points in rows and training curves in
#'   columns.
#' @seealso [farfima()]; [`plot()`][plot.farfima] to draw the curves.
#' @examples
#' sim = sim_farfima(n = 200, d = 0.3, p = 1, q = 0, seed = 1)
#' fit = farfima(sim, d = 0.3, order = c(1, 0), trim = 50, verbose = FALSE)
#' dim(fitted(fit))
#' range(fitted(fit) + residuals(fit) - sim[, 51:200], na.rm = TRUE)
#' @name fitted.farfima
#' @export
fitted.farfima = function(object, scale = c("original", "differenced"), ...) {
  scale = match.arg(scale)
  if (attr(object, "settings")$method == "far") {
    stop("Fitted curves are only available for the FARMA method.", call. = FALSE)
  }
  if (scale == "original") object$fitted else object$fitted_differenced
}
