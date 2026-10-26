# Calibrate the bands from rolling forecasts over the validation sample
calibrate_bands = function(X, X_eval, model, labels, train_end, valid_end, h_max, alpha_set,
                           eta_grid, zeta_grid, verbose, band = c("SR", "NS", "AC"), symmetric = TRUE,
                           integrate = TRUE, KI = NULL, regression = TRUE, scorecast = "none") {
  n_grid = nrow(X_eval)
  E_valid = lapply(seq_len(h_max), function(hh) matrix(numeric(0), nrow = n_grid, ncol = 0))
  F_valid = E_valid
  Y_valid = E_valid
  V_origin = lapply(seq_len(h_max), function(hh) integer(0))
  failures = character(0)

  if (verbose) {
    message(sprintf("Calibrating the %s bands on %d validation origins (%d at h = %d)",
                    paste(band, collapse = ", "), valid_end - train_end, valid_end - train_end - h_max + 1, h_max))
  }
  bar = progress_bar(verbose, valid_end - train_end)

  # Validation errors, forecasts and observed curves for each horizon
  for (t0 in train_end:(valid_end - 1)) {
    if (!is.null(bar)) utils::setTxtProgressBar(bar, t0 - train_end + 1)
    h_here = min(h_max, valid_end - t0)
    pred = origin_forecast(X, t0, h_here, model)

    if (is.character(pred)) {
      failures[labels[t0]] = pred
      next
    }

    truth = X_eval[, (t0 + 1):(t0 + h_here), drop = FALSE]
    for (hh in seq_len(h_here)) {
      E_valid[[hh]] = cbind(E_valid[[hh]], truth[, hh] - pred[, hh])
      F_valid[[hh]] = cbind(F_valid[[hh]], pred[, hh])
      Y_valid[[hh]] = cbind(Y_valid[[hh]], truth[, hh])
      V_origin[[hh]] = c(V_origin[[hh]], t0)
    }
  }
  if (!is.null(bar)) close(bar)

  # SR, NS and AC calibration for each level and horizon
  split_origin = train_end + floor((valid_end - train_end) / 2)
  calibration = list()

  for (alpha in alpha_set) {
    cal = list(gamma = vector("list", h_max), zeta = rep(NA_real_, h_max), Q_hat = rep(NA_real_, h_max),
               scores = vector("list", h_max), eta = 0, n_origin = lengths(V_origin))

    for (hh in seq_len(h_max)) {
      if (ncol(E_valid[[hh]]) == 0) {
        next
      }
      sr = calibrate_zeta_ecp(E_valid = E_valid[[hh]], alpha = alpha, zeta_grid = zeta_grid)
      cal$gamma[[hh]] = sr$gamma
      cal$zeta[hh] = sr$zeta
      cal$scores[[hh]] = conformal_scores(truth = Y_valid[[hh]], forecast = F_valid[[hh]], gamma = sr$gamma)
      cal$Q_hat[hh] = conformal_quantile(cal$scores[[hh]], alpha)
    }

    cal$eta = aci_select_eta(F_by_h = F_valid, Y_by_h = Y_valid, origins_by_h = V_origin,
                             split_origin = split_origin, alpha = alpha, eta_grid = eta_grid)$eta
    cal$aci_err = aci_tuning_errors(F_valid, Y_valid, V_origin, split_origin, alpha, cal$eta)

    # Learning rates of the online bands
    cal$online = list()
    for (m in online_methods(band, symmetric)) {
      sel = online_select_rate(m, symmetric && m != "AcMCP", E_valid, F_valid, Y_valid, V_origin, split_origin, alpha,
                               eta_grid, integrate, KI, regression, scorecast)
      cal$online[[m]] = list(rate = sel$rate, criterion = sel$criterion)
    }
    calibration[[as.character(alpha)]] = cal
  }

  # Validation scores for the online bands
  reference = online_reference(E_valid, V_origin, calibration[[1]]$gamma)

  list(calibration = calibration, failures = failures, reference = reference)
}

# AC error rates on the second half of the validation sample
aci_tuning_errors = function(F_by_h, Y_by_h, origins_by_h, split_origin, alpha, eta) {
  lapply(seq_along(F_by_h), function(hh) {
    o = origins_by_h[[hh]]
    cal = o < split_origin
    tun = o >= split_origin

    if (length(o) == 0 || sum(cal) < 2 || sum(tun) < 1) {
      return(numeric(0))
    }

    Fh = F_by_h[[hh]]
    Yh = Y_by_h[[hh]]
    gamma = aci_gamma(Yh[, cal, drop = FALSE] - Fh[, cal, drop = FALSE])
    scores = conformal_scores(Yh[, cal, drop = FALSE], Fh[, cal, drop = FALSE], gamma)
    o_tun = o[tun]
    idx = which(tun)
    err = rep(NA_real_, length(idx))
    used = rep(FALSE, length(idx))
    a_t = alpha

    for (i in seq_along(idx)) {
      for (j in which(!used & o_tun + hh <= o_tun[i])) {
        if (is.finite(err[j])) {
          a_t = a_t + eta * (alpha - err[j])
        }
        used[j] = TRUE
      }
      band = conformal_interval(Fh[, idx[i]], gamma, aci_quantile(scores, a_t))
      y = Yh[, idx[i]]
      keep = is.finite(y) & is.finite(band$lower) & is.finite(band$upper)
      if (any(keep)) {
        err[i] = mean(y[keep] < band$lower[keep] | y[keep] > band$upper[keep])
      }
    }

    err
  })
}
