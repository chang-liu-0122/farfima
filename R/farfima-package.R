#' farfima: Order Selection and Forecasting in Long-Memory Functional Time Series
#'
#' Fits FARFIMA models to long-memory functional time series, selects their
#' orders and forecasts the curves with conformal prediction bands. Curves are
#' stored as a matrix with grid points in rows and time in columns. `fts`
#' objects from rainbow and ftsa also work.
#'
#' @section Workflow:
#' 1. [split_curves()] splits the curves into training, validation and test
#'    samples, and [plot_curves()] draws them.
#' 2. [farfima()] fits the model. It estimates the memory parameter with
#'    [estimate_memory()], selects the FARMA orders of the fractionally
#'    differenced curves as in [select_order()], estimates their operators as
#'    in [estimate_operators()], and calibrates the conformal prediction bands
#'    on the validation sample. [`fitted()`][fitted.farfima] and
#'    [`residuals()`][fitted.farfima] give the in-sample fit.
#' 3. [`predict()`][predict.farfima] forecasts from the fit, with conformal
#'    prediction bands, and evaluates the forecasts on a test sample.
#' 4. [`plot()`][plot.farfima_forecast] visualises the forecasts, and
#'    [`as.data.frame()`][as.data.frame.farfima_forecast] puts them in a
#'    table.
#'
#' @section Other functions:
#' * [estimate_memory()], [select_order()] and [estimate_operators()]: the
#'   steps of [farfima()], which can also be run on their own.
#'   [select_order()] and [estimate_operators()] take short-memory curves,
#'   such as fractionally differenced curves.
#' * [frac_diff()]: fractional differencing.
#' * [sim_farfima()]: simulated long-memory functional time series.
#'
#' @author Chang Liu, Han Lin Shang, Garth Tarr and Samuel \enc{Müller}{Mueller}.
#'   Maintainer: Chang Liu \email{chang.liu45@students.mq.edu.au}
#' @examples
#' sim = sim_farfima(n = 200, d = 0.3, p = 1, q = 1, no_grid = 51, seed = 1)
#' blocks = split_curves(sim, train = 0.6, valid = 0.2)
#' fit = farfima(blocks$train, blocks$valid, d = "lpwn2", order = "ic_full",
#'                p = 0:2, q = 0:2, h = 5, trim = 30)
#' fc = predict(fit, test = blocks$test, band = "AC")
#' fc
#' plot(fc, h = 1)
#' plot(fc, h = 3, target_test_curve = 10)
#' @importFrom stats ar cov lm.fit lsfit optim optimize predict sd var median coef fitted residuals
"_PACKAGE"
