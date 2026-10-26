# Centre, fractionally difference and trim the curves available at forecast origin t0
origin_history = function(X, t0, d, trim, window, window_size) {
  X_history = X[, 1:t0, drop = FALSE]

  if (window == "expanding") {
    X_mean_sample = X_history
  } else {
    X_mean_sample = X[, max(1, t0 - window_size + 1):t0, drop = FALSE]
  }

  mean_curve = rowMeans(X_mean_sample)
  X_history_centered = sweep(X_history, 1, mean_curve, "-")

  if (d == 0) {
    Y_history = X_history_centered
    Y_estimation = Y_history
  } else {
    Y_history = frac_difference_history(X_history_centered, d)
    if (ncol(Y_history) <= trim) {
      stop("The history is too short after trimming `trim` fractionally differenced curves.", call. = FALSE)
    }
    Y_estimation = Y_history[, (trim + 1):ncol(Y_history), drop = FALSE]
  }

  if (window == "expanding") {
    Y_train = Y_estimation
  } else {
    Y_train = Y_estimation[, max(1, ncol(Y_estimation) - window_size + 1):ncol(Y_estimation), drop = FALSE]
  }

  list(Y_history = Y_history, Y_train = Y_train, mean_curve = mean_curve)
}

# Forecasts from origin t0; check = FALSE keeps an unstable fit
forecast_at_origin = function(X, t0, h, d, p, q, method, Npc, far_pve, trim, window, window_size,
                              check = TRUE) {
  history = origin_history(X, t0, d, trim, window, window_size)
  limit = if (check) c(1, 1, 50) else c(Inf, Inf, Inf)

  out = forecast_original_scale(
    Y_train = history$Y_train,
    d_hat = d,
    h = h,
    method = method,
    p_est = p,
    q_est = q,
    Npc = Npc,
    far_pve = far_pve,
    Y_history_for_integration = history$Y_history,
    farma_ar_radius_limit = limit[1],
    farma_ma_radius_limit = limit[2],
    farma_residual_ratio_limit = limit[3]
  )

  out[[method]]$X_future + matrix(history$mean_curve, nrow = nrow(X), ncol = h)
}

# Forecasts from origin t0, or the error message when the fit fails
origin_forecast = function(X, t0, h, model, check = TRUE) {
  tryCatch(
    forecast_at_origin(X, t0, h, model$d, model$p, model$q, model$method, model$Npc, model$far_pve,
                       model$trim, model$window, model$train_end, check = check),
    error = function(e) conditionMessage(e)
  )
}

# Warn about the forecast origins that failed and were skipped
warn_failures = function(failures, n_origins, sample_name) {
  if (length(failures) == 0) {
    return(invisible(NULL))
  }
  shown = paste(utils::head(names(failures), 6), collapse = ", ")
  if (length(failures) > 6) {
    shown = paste0(shown, " and ", length(failures) - 6, " more")
  }
  warning(sprintf("%d of %d %s forecasts failed and were skipped (origins %s).",
                  length(failures), n_origins, sample_name, shown), call. = FALSE)
}
