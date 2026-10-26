#' Fit a FARFIMA model
#'
#' Fits a functional autoregressive fractionally integrated moving-average
#' (FARFIMA) model to the training sample and calibrates conformal prediction
#' bands on the validation sample.
#'
#' The fit has three steps. The memory parameter is estimated by
#' [estimate_memory()]. The training curves are then centred and fractionally
#' differenced, and the first `trim` differenced curves are dropped, which
#' gives short-memory curves. Their FARMA orders are selected as in
#' [select_order()], and their operators are estimated by the functional
#' Hannan-Rissanen method, as in [estimate_operators()].
#'
#' With a validation sample, forecasts up to `h` steps ahead are made from
#' every validation origin, and their errors calibrate the conformal prediction
#' bands for each horizon and level. A validation origin whose fit fails is
#' skipped with a warning. A warning is also given when the fitted operators are
#' unstable.
#'
#' @param train Training sample: a matrix of curves (grid points by time) or an
#'   `fts` object.
#' @param valid Validation sample, following `train`. It is needed for the
#'   conformal prediction bands and for `order = "cv"`.
#' @param d Memory parameter estimator for [estimate_memory()] (`"lw"`, `"dfa"`,
#'   `"lpwn1"`, `"lpwn2"`, `"bc_lpwn1"` or `"bc_lpwn2"`), or a known value.
#' @param order Order selection method for [select_order()] (`"ic_full"`,
#'   `"ic_auto"`, `"bj"` or `"cv"`), or known orders `c(p, q)`. With `"cv"`,
#'   the forecasts are compared with the validation curves on their original
#'   scale.
#' @param p,q Candidate AR and MA orders searched by the selection method.
#'   Not used when `order` gives the orders. (0, 0) is excluded.
#' @param method `"farma"` for the functional Hannan-Rissanen estimator, or
#'   `"far"` for the FAR(1) benchmark fitted by [far::far()]. The benchmark
#'   ignores `order`, and with `d = 1` it is the FARFIMA(1, 1, 0) benchmark.
#' @param h Maximum forecast horizon of the conformal prediction bands and of
#'   `order = "cv"`. `valid` needs at least `h + 1` curves.
#' @param level Nominal coverage levels of the conformal prediction bands.
#' @param band Conformal prediction bands to calibrate: `"SR"` (scaled
#'   residual), `"NS"` (non-conformity score), `"AC"` (adaptive conformal),
#'   `"PID"` (conformal PID control, proportional-integral-derivative) and
#'   `"AcMCP"` (autocorrelated multistep conformal). See
#'   [`predict()`][predict.farfima].
#' @param window `"expanding"` or `"fixed"` window at each forecast origin. The
#'   fixed window has the length of the training sample.
#' @param valid_eval Optional curves to check the validation forecasts against,
#'   in place of `valid`. The forecasts are still made from `valid`. Use it
#'   when the model is fitted to smoothed curves but the bands should be
#'   calibrated on the raw curves. It must have the same size as `valid`.
#' @param Npc Number of principal components in the operator estimator.
#' @param trim Number of initial fractionally differenced curves dropped before
#'   order selection and estimation. Defaults to \eqn{\min(100, \lfloor T/3
#'   \rfloor)}, where \eqn{T} is the number of training curves.
#' @param far_pve Proportion of variance explained by the components of the
#'   FAR(1) benchmark.
#' @param memory_args List of other arguments to [estimate_memory()], such as
#'   `m_const` or `B`.
#' @param order_args List of other arguments to [select_order()], such as
#'   `n_lags` or `alpha`.
#' @param eta_grid Candidate learning rates of the AC, PID and AcMCP conformal
#'   prediction bands. One rate is chosen for each band and level, by the mean
#'   interval score on the second half of the validation sample.
#' @param zeta_grid Candidate multipliers of the SR band.
#' @param symmetric If `TRUE`, the AC and PID bands are symmetric about the
#'   point forecast. If `FALSE`, their lower and upper edges are calibrated
#'   separately, each at level \eqn{\alpha/2}. The SR and NS bands are always
#'   symmetric and the AcMCP band is always asymmetric.
#' @param integrate Whether the PID and AcMCP bands include the integral term,
#'   which widens the band after repeated under-coverage and narrows it after
#'   repeated over-coverage. In the extreme, the band becomes infinite or
#'   empty.
#' @param KI Scale of the integral term. Defaults to the largest absolute
#'   standardised validation residual at each horizon.
#' @param scorecast Scorecaster of the PID band: `"none"`, `"naive"` (the
#'   latest observed score curve) or `"theta"` (the theta method: exponential
#'   smoothing of the past scores plus half their linear trend).
#' @param acmcp_regression Whether the centre of the AcMCP band also uses the
#'   regression on the scores at shorter horizons. If `FALSE`, only the
#'   MA(h - 1) forecast is used.
#' @param seed Optional seed for the BC-LPWN bootstrap and the white-noise test
#'   of `order = "bj"`.
#' @param verbose Whether to show progress messages and a progress bar.
#' @return An object of class `farfima`, a list with
#' * `d` and `order`: the memory parameter and `c(p = , q = )`;
#' * `phi`, `theta`: the AR and MA operators (grid by grid by lag);
#' * `fitted`, `fitted_differenced`, `residuals`: the in-sample fit on the
#'   scale of the data and on the fractionally differenced scale, and the
#'   residuals, also given by [`fitted()`][fitted.farfima] and
#'   [`residuals()`][fitted.farfima] (`NULL` for the FAR(1) benchmark);
#' * `memory`, `order_selection`: the [estimate_memory()] and [select_order()]
#'   results, or `NULL` when `d` or the orders are given;
#' * `operators`: the [estimate_operators()] result, or `NULL` for the FAR(1)
#'   benchmark;
#' * `calibration`: a table of the calibrated conformal prediction bands by
#'   `level` and `h`:
#'   `SR_zeta` (SR multiplier), `NS_quantile` (NS quantile), `AC_eta`,
#'   `PID_eta` and `AcMCP_eta` (learning rates, for the bands chosen) and
#'   `n_origin` (number of validation origins). `NULL` without a validation
#'   sample.
#'
#' `print()` and `summary()` show the fit, and `coef()` gives `d`, `phi` and
#' `theta`.
#' @seealso [`predict()`][predict.farfima] to forecast from the fit;
#'   [estimate_memory()], [select_order()] and [estimate_operators()] for the
#'   steps on their own.
#' @references
#' Kuenzer, T. (2024). Estimation of functional ARMA models. *Bernoulli*.
#'
#' Li, D., Robinson, P. M. and Shang, H. L. (2020). Long-range dependent curve
#' time series. *Journal of the American Statistical Association*.
#' @examples
#' sim = sim_farfima(n = 200, d = 0.3, p = 1, q = 0, no_grid = 51, seed = 1)
#' blocks = split_curves(sim)
#' fit = farfima(blocks$train, blocks$valid, d = "lpwn2", order = "ic_full",
#'                p = 0:2, q = 0:2, h = 5, trim = 30)
#' fit
#' @export
farfima = function(train,
                    valid = NULL,
                    d = c("bc_lpwn2", "bc_lpwn1", "lpwn2", "lpwn1", "lw", "dfa"),
                    order = c("ic_full", "ic_auto", "bj", "cv"),
                    p = 0:5,
                    q = 0:5,
                    method = c("farma", "far"),
                    h = 10,
                    level = c(0.8, 0.95),
                    band = c("SR", "NS", "AC"),
                    window = c("expanding", "fixed"),
                    valid_eval = NULL,
                    Npc = 1,
                    trim = NULL,
                    far_pve = 0.95,
                    memory_args = list(),
                    order_args = list(),
                    eta_grid = c(0, 0.005, 0.01, 0.02, 0.05, 0.1),
                    zeta_grid = seq(0.1, 50, length.out = 2496),
                    symmetric = TRUE,
                    integrate = TRUE,
                    KI = NULL,
                    scorecast = c("none", "naive", "theta"),
                    acmcp_regression = TRUE,
                    seed = NULL,
                    verbose = TRUE) {
  method = match.arg(method)
  window = match.arg(window)
  band = match.arg(band, c("SR", "NS", "AC", "PID", "AcMCP"), several.ok = TRUE)
  scorecast = match.arg(scorecast)
  train = as_curves(train, "train")
  h = as.integer(h)

  # Check the inputs
  if (length(h) != 1 || is.na(h) || h < 1) {
    stop("`h` must be a single positive integer.", call. = FALSE)
  }
  if (!is.numeric(level) || any(level <= 0 | level >= 1)) {
    stop("`level` must be between 0 and 1.", call. = FALSE)
  }
  if (!is.null(valid)) {
    valid = as_curves(valid, "valid", nrow(train))
    check_follow(train, valid, "the training sample", "valid")
    if (ncol(valid) < h + 1) {
      stop(sprintf("With h = %d, `valid` needs at least %d curves but has %d.", h, h + 1, ncol(valid)), call. = FALSE)
    }
  }
  if (!is.null(valid_eval)) {
    if (is.null(valid)) {
      stop("`valid_eval` needs `valid`.", call. = FALSE)
    }
    valid_eval = as_curves(valid_eval, "valid_eval", nrow(train))
    if (ncol(valid_eval) != ncol(valid)) {
      stop("`valid_eval` must have the same number of curves as `valid`.", call. = FALSE)
    }
  }
  if (method == "far" && !missing(order)) {
    warning("`order` is ignored for the FAR(1) benchmark.", call. = FALSE)
  }
  check_args(memory_args, estimate_memory, "memory_args", c("x", "method", "seed", "verbose"))
  check_args(order_args, select_order, "order_args",
             c("x", "method", "p", "q", "Npc", "valid", "valid_eval", "h", "window", "seed", "verbose"))
  if (is.null(trim)) {
    trim = min(100, floor(ncol(train) / 3))
  }
  check_trim(ncol(train), trim, if (is.numeric(d)) d else 1)

  if (!is.null(seed)) {
    withr::local_seed(seed)
  }

  # Estimate the memory parameter
  memory = NULL
  if (is.numeric(d)) {
    d_hat = as.numeric(d)
  } else {
    memory = do.call(estimate_memory, c(list(x = train, method = match.arg(d), verbose = verbose), memory_args))
    d_hat = memory$d
    if (!isTRUE(attr(memory, "diagnostics")$admissible)) {
      stop(sprintf("The bias-corrected estimate of d is negative (%.4f). Use another estimator or give d.", d_hat),
           call. = FALSE)
    }
  }
  if (length(d_hat) != 1 || !is.finite(d_hat)) {
    stop("`d` must be an estimator name or a single finite value.", call. = FALSE)
  }

  # Select the orders
  order_selection = NULL
  if (method == "far") {
    p_hat = 1L
    q_hat = 0L
  } else if (is.numeric(order)) {
    if (length(order) != 2 || any(order < 0) || all(order == 0)) {
      stop("Known orders must be given as c(p, q) with p + q > 0.", call. = FALSE)
    }
    p_hat = as.integer(order[1])
    q_hat = as.integer(order[2])
  } else {
    order_selection = do.call(order_select, c(
      list(x = train, d = d_hat, trim = trim, method = match.arg(order), p = p, q = q, Npc = Npc,
           valid = valid, valid_eval = valid_eval, h = h, window = window, verbose = verbose),
      order_args
    ))
    p_hat = order_selection$p
    q_hat = order_selection$q
  }

  # Estimate the operators on the trimmed, fractionally differenced training curves
  Y_train = origin_history(train, ncol(train), d_hat, trim, "expanding", ncol(train))$Y_train
  operators = NULL
  fitted_curves = NULL
  if (method == "farma") {
    cols = (if (d_hat == 0) 0 else trim) + seq_len(ncol(Y_train))
    dimnames(Y_train) = list(rownames(train), colnames(train)[cols])
    operators = operator_fit(Y_train, p_hat, q_hat, Npc)
    est = attr(operators, "estimator")
    fitted_curves = train[, cols, drop = FALSE] - operators$residuals
    if (length(attr(operators, "stability")$reasons) > 0) {
      warning(sprintf("The fitted FARFIMA(%d, %.4g, %d) is unstable; try smaller p and q.", p_hat, d_hat, q_hat),
              call. = FALSE)
    }
  } else {
    est = forecast_far_path(Y_train, h = 1, far_pve = far_pve)
  }

  # Calibrate the bands on the validation sample
  calibration = NULL
  if (!is.null(valid)) {
    X = cbind(train, valid)
    X_eval = cbind(train, if (is.null(valid_eval)) valid else valid_eval)
    model = list(d = d_hat, p = p_hat, q = q_hat, method = method, Npc = Npc, far_pve = far_pve,
                 trim = trim, window = window, train_end = ncol(train))
    cal = calibrate_bands(X, X_eval, model, origin_labels(X), ncol(train), ncol(X), h,
                          round(1 - level, 10), eta_grid, zeta_grid, verbose, band, symmetric,
                          integrate, KI, acmcp_regression, scorecast)
    calibration = do.call(rbind, lapply(names(cal$calibration), function(a) {
      cl = cal$calibration[[a]]
      tab = data.frame(level = 1 - as.numeric(a), h = seq_len(h), SR_zeta = cl$zeta, NS_quantile = cl$Q_hat,
                       AC_eta = if (symmetric || is.null(cl$online$AC)) cl$eta else cl$online$AC$rate)
      for (m in intersect(c("PID", "AcMCP"), band)) {
        tab[[paste0(m, "_eta")]] = cl$online[[m]]$rate
      }
      tab$n_origin = cl$n_origin
      tab
    }))
    warn_failures(cal$failures, ncol(valid), "validation")
    n_origin = cal$calibration[[1]]$n_origin
    if (any(n_origin == 0)) {
      stop(sprintf("No validation forecast succeeded at h = %s. Try smaller p and q, or a smaller h.",
                   paste(which(n_origin == 0), collapse = ", ")), call. = FALSE)
    }
  }

  center = rowMeans(train)

  structure(
    list(
      d = d_hat,
      order = c(p = p_hat, q = q_hat),
      phi = operators$phi,
      theta = operators$theta,
      fitted = fitted_curves,
      fitted_differenced = operators$fitted,
      residuals = operators$residuals,
      memory = memory,
      order_selection = order_selection,
      operators = operators,
      calibration = calibration
    ),
    class = "farfima",
    settings = list(method = method, h = h, level = level, band = band, window = window,
                    trim = trim, Npc = Npc, far_pve = far_pve, symmetric = symmetric,
                    integrate = integrate, KI = KI, regression = acmcp_regression, scorecast = scorecast),
    data = list(train = train, valid = valid, valid_eval = valid_eval),
    calibration_detail = if (is.null(valid)) NULL else cal$calibration,
    online_reference = if (is.null(valid)) NULL else cal$reference,
    center = center,
    estimator = est,
    stability = attr(operators, "stability"),
    failures = if (is.null(valid)) character(0) else cal$failures
  )
}
