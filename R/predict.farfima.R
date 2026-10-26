#' Forecast from a FARFIMA fit
#'
#' Point forecasts and conformal prediction bands up to `h` steps beyond the
#' data. With a test sample, forecasts are also made from every test origin and
#' evaluated.
#'
#' The conformal prediction bands are those chosen in [farfima()]:
#' * `"SR"`: scaled residual band, the forecast plus and minus a multiplier of
#'   the point-wise scale function;
#' * `"NS"`: non-conformity score band, with the width set by a quantile of the
#'   standardised non-conformity scores;
#' * `"AC"`: adaptive conformal band, an NS band whose quantile level is
#'   updated over the test sample;
#' * `"PID"`: conformal PID control band, whose multiplier is updated by a
#'   proportional term, an integral term and an optional scorecaster;
#' * `"AcMCP"`: autocorrelated multistep conformal band, an asymmetric PID band
#'   shifted by a forecast of the standardised residual curve.
#'
#' The AC, PID and AcMCP bands are updated from the error rates of earlier
#' bands, once their targets are observed.
#'
#' In the test evaluation, the memory parameter and the orders stay fixed, and
#' the operators and the mean are re-estimated at each origin. An origin whose
#' fit fails is skipped with a warning. The forecast beyond the data starts
#' after the last test curve, or after the last validation curve when there is
#' no test sample.
#'
#' @param object A fit from [farfima()].
#' @param test Optional test sample, following the validation sample.
#' @param h Maximum forecast horizon, at most the `h` of the fit.
#' @param band Conformal prediction bands to return, from those calibrated in
#'   [farfima()].
#' @param test_eval Optional curves to check the test forecasts against, in
#'   place of `test`. The forecasts are still made from `test`.
#' @param verbose Whether to show progress messages and a progress bar.
#' @param ... Not used.
#' @return An object of class `farfima_forecast`, a list with
#' * `point_forecast`: the forecast curves beyond the data (grid by `h`);
#' * `bands`: the conformal prediction bands, with the `lower` and `upper`
#'   matrices (grid by `h`) for each band and level (such as `"80%"`);
#' * `evaluation`, with a test sample: the tables `point` (RMSPE, the root mean
#'   squared prediction error, and MAPE, the mean absolute prediction error),
#'   `interval` (ECP, CPD and MIS of the bands) and `wis` (weighted interval
#'   score).
#' @seealso [`plot()`][plot.farfima_forecast] and
#'   [`as.data.frame()`][as.data.frame.farfima_forecast] for the results.
#' @references
#' Gibbs, I. and Candes, E. J. (2021). Adaptive conformal inference under
#' distribution shift. *Advances in Neural Information Processing Systems*.
#'
#' Angelopoulos, A. N., Candes, E. J. and Tibshirani, R. J. (2023). Conformal
#' PID control for time series prediction. *Advances in Neural Information
#' Processing Systems*.
#'
#' Wang, X. and Hyndman, R. J. (2024). Online conformal inference for
#' multi-step time series forecasting. arXiv:2410.13115.
#'
#' Wang, X. and Girolimetto, D. (2026). conformalForecast: an R package for
#' sequential conformal prediction in time series. arXiv:2610.21648.
#'
#' Bracher, J., Ray, E. L., Gneiting, T. and Reich, N. G. (2021). Evaluating
#' epidemic forecasts in an interval format. *PLOS Computational Biology*.
#' @examples
#' sim = sim_farfima(n = 200, d = 0.3, p = 1, q = 0, no_grid = 51, seed = 1)
#' blocks = split_curves(sim)
#' fit = farfima(blocks$train, blocks$valid, d = 0.3, order = c(1, 0), h = 5, trim = 30)
#' fc = predict(fit, test = blocks$test)
#' fc
#' @export
predict.farfima = function(object,
                            test = NULL,
                            h = attr(object, "settings")$h,
                            band = attr(object, "settings")$band,
                            test_eval = NULL,
                            verbose = TRUE,
                            ...) {
  moved = intersect(names(list(...)), c("valid", "level", "window", "valid_eval", "eta_grid", "zeta_grid"))
  if (length(moved) > 0) {
    stop(sprintf("`%s` is set in farfima(), not predict().", moved[1]), call. = FALSE)
  }

  # Check the inputs
  h_max = as.integer(h)
  set = attr(object, "settings")
  dat = attr(object, "data")
  calibrated = !is.null(attr(object, "calibration_detail"))
  band = match.arg(band, c("SR", "NS", "AC", "PID", "AcMCP"), several.ok = TRUE)
  if (any(!band %in% set$band)) {
    stop(sprintf("Band %s was not chosen in farfima().", paste(setdiff(band, set$band), collapse = ", ")),
         call. = FALSE)
  }
  if (length(h_max) != 1 || is.na(h_max) || h_max < 1) {
    stop("`h` must be a single positive integer.", call. = FALSE)
  }
  if (calibrated && h_max > set$h) {
    stop(sprintf("The bands were calibrated for h <= %d in farfima().", set$h), call. = FALSE)
  }

  train = dat$train
  valid = dat$valid
  if (!is.null(test)) {
    test = as_curves(test, "test", nrow(train))
    if (is.null(valid)) {
      check_follow(train, test, "the training sample", "test")
    } else {
      check_follow(valid, test, "the validation sample", "test")
    }
  }
  if (!is.null(test_eval)) {
    if (is.null(test)) {
      stop("`test_eval` needs `test`.", call. = FALSE)
    }
    test_eval = as_curves(test_eval, "test_eval", nrow(train))
    if (ncol(test_eval) != ncol(test)) {
      stop("`test_eval` must have the same number of curves as `test`.", call. = FALSE)
    }
  }

  X = cbind(train, valid, test)
  X_eval = cbind(train,
                 if (is.null(dat$valid_eval)) valid else dat$valid_eval,
                 if (is.null(test_eval)) test else test_eval)
  alpha_set = round(1 - set$level, 10)
  T_all = ncol(X)
  n_grid = nrow(X)
  train_end = ncol(train)
  valid_end = train_end + if (is.null(valid)) 0L else ncol(valid)
  has_time = !is.null(colnames(X)) && all(nzchar(colnames(X)))
  labels = origin_labels(X)
  model = list(d = object$d, p = object$order[["p"]], q = object$order[["q"]], method = set$method,
               Npc = set$Npc, far_pve = set$far_pve, trim = set$trim, window = set$window,
               train_end = train_end)
  calibration = attr(object, "calibration_detail")
  aci_err = if (calibrated) lapply(calibration, function(cal) cal$aci_err) else NULL
  evaluation = NULL
  test_forecasts = NULL
  failures = character(0)
  symmetric = !isFALSE(set$symmetric)
  online = online_methods(band, symmetric)
  F_run = lapply(seq_len(h_max), function(hh) matrix(numeric(0), nrow = n_grid, ncol = 0))
  Y_run = F_run
  O_run = lapply(seq_len(h_max), function(hh) integer(0))
  rec_index = integer(0)

  # Point forecasts beyond the data
  mean_future = origin_forecast(X, T_all, h_max, model)
  if (is.character(mean_future) && grepl("unstable FARMA fit", mean_future)) {
    warning("The forecast beyond the data comes from an unstable FARMA fit.", call. = FALSE)
    mean_future = origin_forecast(X, T_all, h_max, model, check = FALSE)
  }
  if (is.character(mean_future)) {
    stop(sprintf("The forecast beyond the data failed: %s", mean_future), call. = FALSE)
  }
  dimnames(mean_future) = list(rownames(train), future_labels(colnames(X), h_max))

  # Rolling forecasts over the test sample
  if (!is.null(test)) {
    test_steps = valid_end:(T_all - 1)
    sse = numeric(h_max)
    sae = numeric(h_max)
    n_pt = numeric(h_max)
    test_forecasts = list()
    acc = list()
    state = list()

    for (alpha in alpha_set) {
      a_name = as.character(alpha)
      acc[[a_name]] = lapply(seq_len(h_max), function(hh) {
        list(SR = c(0, 0, 0), NS = c(0, 0, 0), AC = c(0, 0, 0), PID = c(0, 0, 0), AcMCP = c(0, 0, 0))
      })
      state[[a_name]] = lapply(seq_len(h_max), function(hh) {
        list(a_t = alpha, origin = integer(0), err = numeric(0), used = logical(0))
      })
    }

    if (verbose) {
      message(sprintf("Evaluating on %d test origins", length(test_steps)))
    }
    bar = progress_bar(verbose, length(test_steps))

    for (t0 in test_steps) {
      if (!is.null(bar)) utils::setTxtProgressBar(bar, t0 - valid_end + 1)
      h_here = min(h_max, T_all - t0)
      pred = origin_forecast(X, t0, h_here, model)
      if (is.character(pred)) {
        failures[labels[t0]] = pred
        next
      }

      # Point forecast errors
      truth = X_eval[, (t0 + 1):(t0 + h_here), drop = FALSE]
      err = pred - truth
      ok = is.finite(err)
      err[!ok] = 0
      sse[seq_len(h_here)] = sse[seq_len(h_here)] + colSums(err^2)
      sae[seq_len(h_here)] = sae[seq_len(h_here)] + colSums(abs(err))
      n_pt[seq_len(h_here)] = n_pt[seq_len(h_here)] + colSums(ok)
      record = list(origin = labels[t0], target = labels[(t0 + 1):(t0 + h_here)],
                    target_test_curve = (t0 + 1):(t0 + h_here) - valid_end,
                    point_forecast = pred, observed = truth, bands = list())
      for (hh in seq_len(h_here)) {
        F_run[[hh]] = cbind(F_run[[hh]], pred[, hh])
        Y_run[[hh]] = cbind(Y_run[[hh]], truth[, hh])
        O_run[[hh]] = c(O_run[[hh]], t0)
      }
      rec_index[as.character(t0)] = length(test_forecasts) + 1

      # SR, NS and AC bands at this origin
      for (alpha in alpha_set) {
        if (!calibrated) {
          next
        }
        a_name = as.character(alpha)
        cal = calibration[[a_name]]
        lo = list(SR = pred * NA, NS = pred * NA, AC = pred * NA)
        up = lo

        for (hh in seq_len(h_here)) {
          gamma_h = as.numeric(cal$gamma[[hh]])
          f_h = pred[, hh]
          y_h = truth[, hh]
          if (length(gamma_h) != length(f_h)) {
            next
          }
          a = acc[[a_name]][[hh]]

          if (is.finite(cal$zeta[hh])) {
            hw = cal$zeta[hh] * gamma_h
            e_h = f_h - y_h
            keep = is.finite(e_h) & is.finite(hw)
            if (any(keep)) {
              a$SR = a$SR + c(sum(abs(e_h[keep]) <= hw[keep]),
                              sum(2 * hw[keep] + 2 / alpha * pmax(abs(e_h[keep]) - hw[keep], 0)),
                              sum(keep))
            }
            lo$SR[, hh] = f_h - hw
            up$SR[, hh] = f_h + hw
          }

          band_ns = conformal_interval(f_h, gamma_h, cal$Q_hat[hh])
          a$NS = a$NS + band_tally(y_h, band_ns, alpha)
          lo$NS[, hh] = band_ns$lower
          up$NS[, hh] = band_ns$upper

          if (!symmetric) {
            acc[[a_name]][[hh]] = a
            next
          }
          st = state[[a_name]][[hh]]
          for (j in which(!st$used & st$origin + hh <= t0)) {
            st$a_t = st$a_t + cal$eta * (alpha - st$err[j])
            st$used[j] = TRUE
          }
          band_ac = conformal_interval(f_h, gamma_h, aci_quantile(cal$scores[[hh]], st$a_t))
          a$AC = a$AC + band_tally(y_h, band_ac, alpha)
          lo$AC[, hh] = band_ac$lower
          up$AC[, hh] = band_ac$upper
          keep = is.finite(y_h) & is.finite(band_ac$lower) & is.finite(band_ac$upper)
          if (any(keep)) {
            st$origin = c(st$origin, t0)
            st$err = c(st$err, mean(y_h[keep] < band_ac$lower[keep] | y_h[keep] > band_ac$upper[keep]))
            st$used = c(st$used, FALSE)
          }

          state[[a_name]][[hh]] = st
          acc[[a_name]][[hh]] = a
        }

        for (m in band) {
          if (m %in% online) {
            record$bands[[m]][[level_label(1 - alpha)]] = list(lower = pred * NA, upper = pred * NA)
          } else {
            record$bands[[m]][[level_label(1 - alpha)]] = list(lower = lo[[m]], upper = up[[m]])
          }
        }
      }

      test_forecasts[[length(test_forecasts) + 1]] = record
    }
    if (!is.null(bar)) close(bar)
    warn_failures(failures, length(test_steps), "test")

    if (calibrated) {
      aci_err = lapply(state, function(s) lapply(s, function(st) st$err))
    }
  }

  # Online bands over the test origins
  online_runs = list()
  if (calibrated && length(online) > 0) {
    F_m = lapply(seq_len(h_max), function(hh) cbind(F_run[[hh]], mean_future[, hh]))
    Y_m = lapply(seq_len(h_max), function(hh) cbind(Y_run[[hh]], NA))
    O_m = lapply(seq_len(h_max), function(hh) c(O_run[[hh]], T_all))
    for (alpha in alpha_set) {
      a_name = as.character(alpha)
      for (m in online) {
        online_runs[[a_name]][[m]] = online_bands(m, symmetric && m != "AcMCP", F_m, Y_m, O_m,
                                                  attr(object, "online_reference"), alpha,
                                                  calibration[[a_name]]$online[[m]]$rate,
                                                  !isFALSE(set$integrate), set$KI, !isFALSE(set$regression),
                                                  if (is.null(set$scorecast)) "none" else set$scorecast)
      }
    }
  }

  if (!is.null(test)) {
    for (a_name in names(online_runs)) {
      lev = level_label(1 - as.numeric(a_name))
      for (m in names(online_runs[[a_name]])) {
        for (hh in seq_len(h_max)) {
          run = online_runs[[a_name]][[m]][[hh]]
          if (is.null(run)) {
            next
          }
          acc[[a_name]][[hh]][[m]] = run$tally
          for (i in seq_len(length(run$origins) - 1)) {
            k = rec_index[[as.character(run$origins[i])]]
            test_forecasts[[k]]$bands[[m]][[lev]]$lower[, hh] = run$lower[, i]
            test_forecasts[[k]]$bands[[m]][[lev]]$upper[, hh] = run$upper[, i]
          }
        }
      }
    }

    # Evaluation tables
    point = data.frame(h = seq_len(h_max), RMSPE = sqrt(sse / n_pt), MAPE = sae / n_pt)
    interval = NULL
    wis = NULL
    if (calibrated) {
      interval = do.call(rbind, lapply(alpha_set, function(alpha) {
        do.call(rbind, lapply(band, function(m) {
          do.call(rbind, lapply(seq_len(h_max), function(hh) {
            v = acc[[as.character(alpha)]][[hh]][[m]]
            ecp = if (v[3] > 0) v[1] / v[3] else NA_real_
            data.frame(band = m, level = 1 - alpha, h = hh, ECP = ecp,
                       CPD = abs(ecp - (1 - alpha)), MIS = if (v[3] > 0) v[2] / v[3] else NA_real_)
          }))
        }))
      }))

      wis = do.call(rbind, lapply(band, function(m) {
        do.call(rbind, lapply(seq_len(h_max), function(hh) {
          rows = interval[interval$band == m & interval$h == hh, ]
          data.frame(band = m, h = hh,
                     WIS = (0.5 * point$MAPE[hh] + sum((1 - rows$level) / 2 * rows$MIS)) / (nrow(rows) + 0.5))
        }))
      }))
    }
    evaluation = list(point = point, interval = interval, wis = wis)
  }

  if (verbose && !calibrated) {
    message("No validation sample in the fit: point forecasts only.")
  }

  # Prediction bands of the forecast beyond the data
  bands = list()
  if (calibrated) {
    for (alpha in alpha_set) {
      a_name = as.character(alpha)
      cal = calibration[[a_name]]
      lev = level_label(1 - alpha)
      lower = upper = list(SR = mean_future * NA, NS = mean_future * NA, AC = mean_future * NA,
                           PID = mean_future * NA, AcMCP = mean_future * NA)

      for (hh in seq_len(h_max)) {
        gamma_h = as.numeric(cal$gamma[[hh]])
        if (length(gamma_h) != n_grid) {
          next
        }
        f_h = mean_future[, hh]
        lower$SR[, hh] = f_h - cal$zeta[hh] * gamma_h
        upper$SR[, hh] = f_h + cal$zeta[hh] * gamma_h
        b = conformal_interval(f_h, gamma_h, cal$Q_hat[hh])
        lower$NS[, hh] = b$lower
        upper$NS[, hh] = b$upper
        if (symmetric) {
          err_h = aci_err[[a_name]][[hh]]
          a_t = alpha + cal$eta * sum(alpha - err_h[is.finite(err_h)])
          b = conformal_interval(f_h, gamma_h, aci_quantile(cal$scores[[hh]], a_t))
          lower$AC[, hh] = b$lower
          upper$AC[, hh] = b$upper
        }
        for (m in names(online_runs[[a_name]])) {
          run = online_runs[[a_name]][[m]][[hh]]
          if (!is.null(run)) {
            lower[[m]][, hh] = run$lower[, ncol(run$lower)]
            upper[[m]][, hh] = run$upper[, ncol(run$upper)]
          }
        }
      }

      for (m in band) {
        bands[[m]][[lev]] = list(lower = lower[[m]], upper = upper[[m]])
      }
    }
  }

  structure(
    list(point_forecast = mean_future, bands = bands, evaluation = evaluation),
    class = "farfima_forecast",
    origin = labels[T_all],
    settings = list(h = h_max, level = set$level, band = if (calibrated) band else character(0),
                    window = set$window),
    model = object,
    curves = X,
    test_forecasts = test_forecasts,
    n_test = if (is.null(test)) 0L else ncol(test),
    last_curve = X[, T_all],
    has_time = has_time,
    failures = failures
  )
}
