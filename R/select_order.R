#' Select the FARMA orders
#'
#' Selects the FARMA(p, q) orders of short-memory curves. Long-memory curves
#' need to be fractionally differenced first, for example with [frac_diff()].
#' [farfima()] does this on the centred training curves and drops the first
#' `trim` before selecting the orders.
#'
#' * `"ic_full"` (functional IC) minimises the adaptive functional
#'   Hannan-Quinn information criterion over all candidates.
#' * `"ic_auto"` (auto functional IC) uses the same adaptive functional
#'   Hannan-Quinn information criterion with a stepwise search from a few
#'   starting models. Its noise level comes from the largest model.
#' * `"bj"` (functional Box-Jenkins) tries the candidates by increasing
#'   `p + q` and stops at the first whose residuals pass the functional
#'   portmanteau test.
#' * `"cv"` (cross-validation) minimises the RMSPE, or MAPE, of the forecasts
#'   from the validation origins, averaged over horizons. A candidate whose fit
#'   fails at any origin is dropped.
#'
#' @param x A matrix of short-memory training curves (grid points by time) or
#'   an `fts` object. For `plot()`, a `farfima_order` object.
#' @param method `"ic_full"`, `"ic_auto"`, `"bj"` or `"cv"`.
#' @param p,q Candidate AR and MA orders. (0, 0) is excluded.
#' @param Npc Number of principal components in the operator estimator.
#' @param n_lags Number of residual autocovariance lags. Defaults to
#'   \eqn{\lfloor \log T \rfloor} for the criteria and 20 for `"bj"`.
#' @param alpha Significance level of the white-noise test (`"bj"`).
#' @param iid Whether the white-noise test assumes iid noise rather than weak
#'   white noise (`"bj"`).
#' @param starts Starting models, a data frame with columns `p` and `q`
#'   (`"ic_auto"`). Defaults to (2, 2), (1, 0) and (0, 1).
#' @param valid Validation curves following `x`, also short-memory (`"cv"`).
#' @param valid_eval Optional curves to check the validation forecasts against,
#'   in place of `valid` (`"cv"`).
#' @param h Maximum forecast horizon (`"cv"`).
#' @param window `"expanding"` or `"fixed"` window (`"cv"`).
#' @param cv_score `"rmse"` or `"mae"` (`"cv"`).
#' @param seed Optional seed. The white-noise test of `"bj"` uses Monte Carlo
#'   integration, so its p-values change between runs without one.
#' @param verbose Whether to show progress messages and a progress bar.
#' @param main Plot title.
#' @param ... Other arguments passed to [graphics::image()], or to
#'   [graphics::barplot()] for `"bj"`.
#' @return An object of class `farfima_order` with the selected `p` and `q`
#'   and a `table` of the candidates. `"ic_auto"` also returns its
#'   `search_path`. `plot()` shows the criterion or CV score of each candidate
#'   with the selected model boxed, or the p-values for `"bj"`.
#' @seealso [estimate_operators()] for the operators, and [farfima()] for the
#'   whole fit.
#' @references
#' Hannan, E. J. and Quinn, B. G. (1979). The determination of the order of an
#' autoregression. *Journal of the Royal Statistical Society: Series B*.
#'
#' Kokoszka, P., Rice, G. and Shang, H. L. (2017). Inference for the
#' autocovariance of a functional time series under conditional
#' heteroscedasticity. *Journal of Multivariate Analysis*.
#'
#' Kuenzer, T. (2024). Estimation of functional ARMA models. *Bernoulli*.
#' @examples
#' x = sim_farfima(n = 300, d = 0, p = 1, q = 0, seed = 1)
#' o = select_order(x, method = "ic_full", p = 0:2, q = 0:2)
#' o
#' plot(o)
#' @export
select_order = function(x,
                        method = c("ic_full", "ic_auto", "bj", "cv"),
                        p = 0:5,
                        q = 0:5,
                        Npc = 1,
                        n_lags = NULL,
                        alpha = 0.05,
                        iid = FALSE,
                        starts = NULL,
                        valid = NULL,
                        valid_eval = NULL,
                        h = 10,
                        window = c("expanding", "fixed"),
                        cv_score = c("rmse", "mae"),
                        seed = NULL,
                        verbose = FALSE) {
  order_select(x, d = 0, trim = 0, method = match.arg(method), p = p, q = q, Npc = Npc, n_lags = n_lags,
               alpha = alpha, iid = iid, starts = starts, valid = valid, valid_eval = valid_eval, h = h,
               window = match.arg(window), cv_score = match.arg(cv_score), seed = seed, verbose = verbose)
}
