# Smallest SR multiplier reaching the nominal coverage
calibrate_zeta_ecp = function(E_valid,
                               alpha,
                               zeta_grid = seq(0.1, 50, length.out = 2496))
{
  E_valid = as.matrix(E_valid)

  E_centered = t(scale(t(E_valid), center = TRUE, scale = FALSE))

  gamma = apply(E_centered, 1, sd, na.rm = TRUE)

  gamma[!is.finite(gamma) | gamma <= 0] = NA_real_

  ecp_grid = sapply(zeta_grid, function(z)
  {
    threshold = matrix(
      z * gamma,
      nrow = nrow(E_valid),
      ncol = ncol(E_valid),
      byrow = FALSE
    )

    mean(abs(E_valid) <= threshold, na.rm = TRUE)
  })

  cand = zeta_grid[ecp_grid >= (1 - alpha)]

  if (length(cand) == 0) {
    zeta = NA_real_
  } else {
    zeta = min(cand)
  }

  list(
    gamma = gamma,
    zeta = zeta,
    ecp_grid = ecp_grid,
    zeta_grid = zeta_grid,
    target = 1 - alpha,
    n_valid = ncol(E_valid),
    n_grid = nrow(E_valid)
  )
}

# Forecast plus and minus the point-wise scale function
conformal_base_interval = function(forecast, gamma) {
  lower = forecast - gamma
  upper = forecast + gamma

  list(lower = lower, upper = upper, width = upper - lower)
}

# Width-standardised non-conformity scores, sorted
conformal_scores = function(truth, forecast, gamma) {
  base = conformal_base_interval(forecast, gamma)

  s = pmax(base$lower - truth, truth - base$upper) / base$width

  sort(s[is.finite(s)])
}

# Calibrated quantile of the non-conformity scores
conformal_quantile = function(scores_sorted, alpha) {
  n = length(scores_sorted)

  if (n == 0L) {
    return(NA_real_)
  }

  k = ceiling((n + 1) * (1 - alpha))

  if (k > n) {
    return(Inf)
  }

  scores_sorted[k]
}

# Non-conformity score band for a given quantile
conformal_interval = function(forecast, gamma, Q) {
  base = conformal_base_interval(forecast, gamma)

  list(lower = base$lower - Q * base$width,
       upper = base$upper + Q * base$width)
}

# Interval score of each observation
interval_score_cells = function(truth, lower, upper, alpha) {
  (upper - lower) +
    (2 / alpha) * (lower - truth) * (truth < lower) +
    (2 / alpha) * (truth - upper) * (truth > upper)
}

# Point-wise scale function from validation errors
aci_gamma = function(E) {
  E = as.matrix(E)

  gamma = apply(t(scale(t(E), center = TRUE, scale = FALSE)), 1, sd, na.rm = TRUE)

  gamma[!is.finite(gamma) | gamma <= 0] = NA_real_

  gamma
}

# Non-conformity score quantile at the adaptive level a_t
aci_quantile = function(scores_sorted, a_t) {
  n = length(scores_sorted)

  if (n == 0L || !is.finite(a_t)) {
    return(NA_real_)
  }

  k = ceiling((n + 1) * (1 - a_t))

  if (k < 1L) {
    return(-0.5)
  }

  scores_sorted[min(k, n)]
}

# Run the adaptive conformal updates over the tuning origins and score the bands
aci_run_sequence = function(scores_sorted, gamma, forecast, truth, origins, h,
                             alpha, eta) {
  n_seq = length(origins)
  err = rep(NA_real_, n_seq)
  used = rep(FALSE, n_seq)
  a_t = alpha

  covered = 0
  n_cells = 0
  score_sum = 0

  for (i in seq_len(n_seq)) {
    for (j in which(!used & origins + h <= origins[i])) {
      if (is.finite(err[j])) {
        a_t = a_t + eta * (alpha - err[j])
      }

      used[j] = TRUE
    }

    band = conformal_interval(forecast[[i]], gamma, aci_quantile(scores_sorted, a_t))
    y = truth[[i]]
    keep = is.finite(y) & is.finite(band$lower) & is.finite(band$upper)

    if (!any(keep)) {
      next
    }

    inside = y[keep] >= band$lower[keep] & y[keep] <= band$upper[keep]
    err[i] = mean(!inside)

    covered = covered + sum(inside)
    n_cells = n_cells + sum(keep)
    score_sum = score_sum +
      sum(interval_score_cells(y[keep], band$lower[keep], band$upper[keep], alpha))
  }

  c(covered = covered, n = n_cells, score_sum = score_sum)
}

# AC learning rate with the smallest mean interval score
aci_select_eta = function(F_by_h, Y_by_h, origins_by_h, split_origin, alpha,
                           eta_grid = c(0, 0.005, 0.01, 0.02, 0.05, 0.1)) {
  eta_grid = sort(unique(eta_grid))

  prepared = list()

  for (h in seq_along(F_by_h)) {
    Fh = F_by_h[[h]]
    Yh = Y_by_h[[h]]
    o = origins_by_h[[h]]

    if (is.null(Fh) || length(o) == 0L || ncol(Fh) != length(o)) {
      next
    }

    cal = o < split_origin
    tun = o >= split_origin

    if (sum(cal) < 2L || sum(tun) < 1L) {
      next
    }

    gamma = aci_gamma(Yh[, cal, drop = FALSE] - Fh[, cal, drop = FALSE])
    scores = conformal_scores(Yh[, cal, drop = FALSE], Fh[, cal, drop = FALSE], gamma)

    if (length(scores) == 0L) {
      next
    }

    prepared[[length(prepared) + 1L]] = list(
      h = h,
      gamma = gamma,
      scores = scores,
      forecast = lapply(which(tun), function(c) Fh[, c]),
      truth = lapply(which(tun), function(c) Yh[, c]),
      origins = o[tun]
    )
  }

  # Mean interval score of each learning rate over the horizons
  criterion = vapply(eta_grid, function(eta) {
    is_h = vapply(prepared, function(x) {
      r = aci_run_sequence(x$scores, x$gamma, x$forecast, x$truth,
                            x$origins, x$h, alpha, eta)
      if (r[["n"]] > 0) r[["score_sum"]] / r[["n"]] else NA_real_
    }, numeric(1))

    if (any(is.finite(is_h))) mean(is_h, na.rm = TRUE) else NA_real_
  }, numeric(1))

  if (!any(is.finite(criterion))) {
    return(list(eta = 0, criterion = criterion, eta_grid = eta_grid))
  }

  rounded = round(criterion, 10)
  best = which(rounded == min(rounded, na.rm = TRUE))[1]

  list(eta = eta_grid[best], criterion = criterion, eta_grid = eta_grid)
}


# Bands updated by online_bands()
online_methods = function(band, symmetric) {
  intersect(band, c(if (!symmetric) "AC", "PID", "AcMCP"))
}

# Quantile of sorted scores, capped at the smallest and largest
capped_quantile = function(sorted, a) {
  n = length(sorted)
  if (n == 0 || !is.finite(a)) {
    return(NA_real_)
  }
  sorted[min(max(ceiling((n + 1) * (1 - a)), 1), n)]
}

# PID integral term KI tan(x log n / (Csat n))
pid_integrator = function(x, n, KI) {
  Csat = if (n > 1) 2 / pi * (ceiling(log(n) * 0.01) - 1 / log(n)) else NA_real_
  if (KI == 0 || !is.finite(Csat) || Csat <= 0) {
    return(0)
  }
  arg = x * log(n) / (Csat * n)
  if (arg >= pi / 2) {
    return(Inf)
  }
  if (arg <= -pi / 2) {
    return(-Inf)
  }
  KI * tan(arg)
}

# Standardised validation residuals and their leading principal component
online_reference = function(E_by_h, origins_by_h, gamma_by_h = NULL) {
  ref = lapply(seq_along(E_by_h), function(hh) {
    E = E_by_h[[hh]]
    if (is.null(E) || ncol(E) < 2) {
      return(NULL)
    }
    gamma = if (is.null(gamma_by_h)) aci_gamma(E) else as.numeric(gamma_by_h[[hh]])
    if (length(gamma) != nrow(E) || any(!is.finite(gamma))) {
      return(NULL)
    }
    list(gamma = gamma, z = E / gamma, origins = origins_by_h[[hh]])
  })
  Z = do.call(cbind, lapply(ref, function(r) r$z))
  basis = if (is.null(Z) || ncol(Z) < 2) NULL else stats::prcomp(t(Z), center = TRUE, scale. = FALSE)$rotation[, 1]
  list(h = ref, basis = basis)
}

# Leading principal component scores of the residual curves, named by origin
pc_scores = function(F_by_h, Y_by_h, origins_by_h, reference) {
  lapply(seq_along(F_by_h), function(hh) {
    rf = reference$h[[hh]]
    if (is.null(rf)) {
      return(numeric(0))
    }
    s = stats::setNames(as.numeric(crossprod(reference$basis, rf$z)), rf$origins)
    if (!is.null(F_by_h[[hh]]) && ncol(F_by_h[[hh]]) > 0) {
      z = (Y_by_h[[hh]] - F_by_h[[hh]]) / rf$gamma
      ok = colSums(!is.finite(z)) == 0
      s = c(s, stats::setNames(as.numeric(crossprod(reference$basis, z[, ok, drop = FALSE])), origins_by_h[[hh]][ok]))
    }
    s
  })
}

# Intercept of an MA(h - 1) fit to the scores
ma_intercept = function(x, h) {
  x = x[is.finite(x)]
  if (h == 1) {
    return(mean(x))
  }
  if (length(x) < h + 2) {
    return(NA_real_)
  }
  fit = tryCatch(suppressWarnings(stats::arima(x, order = c(0, 0, h - 1), method = "CSS-ML")), error = function(e) NULL)
  if (is.null(fit)) NA_real_ else unname(fit$coef[["intercept"]])
}

# Regression of the h-step score on the shorter-horizon scores
cross_horizon_score = function(scores, hh, t, plug) {
  known = Reduce(intersect, lapply(scores[seq_len(hh)], names))
  known = known[as.numeric(known) + hh <= t]
  if (length(known) < hh + 2 || !all(is.finite(plug))) {
    return(NA_real_)
  }
  design = as.data.frame(sapply(seq_len(hh), function(jj) scores[[jj]][known]))
  names(design) = paste0("V", seq_len(hh))
  fit = stats::lm(stats::as.formula(paste0("V", hh, " ~ .")), data = design)
  new = stats::setNames(as.data.frame(as.list(plug)), paste0("V", seq_len(hh - 1)))
  as.numeric(stats::predict(fit, newdata = new))
}

# Naive or theta forecast of the scores at each grid point
score_forecast = function(Z, k, type) {
  n = ncol(Z)
  if (type == "naive" || n < 3) {
    return(Z[, n])
  }
  tt = seq_len(n) - 1
  slope = as.numeric(Z %*% (tt - mean(tt))) / sum((tt - mean(tt))^2)
  best_sse = rep(Inf, nrow(Z))
  best_level = Z[, n]
  best_a = rep(1, nrow(Z))

  # Closed-form initial level for each smoothing weight
  for (a in seq(0.01, 0.99, by = 0.01)) {
    w = (1 - a)^tt
    e0 = matrix(0, nrow(Z), n)
    level = 0
    for (i in seq_len(n)) {
      e0[, i] = Z[, i] - level
      level = level + a * e0[, i]
    }
    l0 = as.numeric(e0 %*% w) / sum(w^2)
    sse = rowSums((e0 - outer(l0, w))^2)
    better = is.finite(sse) & sse < best_sse
    best_sse[better] = sse[better]
    best_level[better] = (level + (1 - a)^n * l0)[better]
    best_a[better] = a
  }
  best_level + slope / 2 * (k - 1 + (1 - (1 - best_a)^n) / best_a)
}

# AC, PID or AcMCP bands over a run of forecast origins
online_bands = function(method, symmetric, F_by_h, Y_by_h, origins_by_h, reference, alpha, rate,
                        integrate = TRUE, KI = NULL, regression = TRUE, scorecast = "none") {
  out = vector("list", length(F_by_h))
  symmetric = symmetric && method != "AcMCP"
  use_scorecast = method == "PID" && scorecast != "none"
  basis = reference$basis
  use_regression = method == "AcMCP" && regression && !is.null(basis)
  scores_pc = if (use_regression) pc_scores(F_by_h, Y_by_h, origins_by_h, reference) else NULL

  for (hh in seq_along(F_by_h)) {
    rf = reference$h[[hh]]
    Fh = F_by_h[[hh]]
    Yh = Y_by_h[[hh]]
    o = origins_by_h[[hh]]
    if (is.null(rf) || is.null(Fh) || ncol(Fh) == 0) {
      next
    }
    gamma = rf$gamma
    z_cal = rf$z

    # One multiplier for a symmetric band, one per edge otherwise
    sides = if (symmetric) "both" else c("lower", "upper")
    target = if (symmetric) alpha else alpha / 2
    track = if (symmetric) abs(z_cal) else z_cal

    # Validation scores less their centre
    r = track
    if (method == "AcMCP") {
      r = track - rowMeans(track)
    } else if (use_scorecast) {
      o_cal = rf$origins
      past = sapply(seq_along(o_cal), function(cc) {
        prev = which(o_cal <= o_cal[cc] - hh)
        if (length(prev) == 0) rep(NA_real_, nrow(track)) else score_forecast(track[, prev, drop = FALSE], hh, scorecast)
      })
      if (any(is.finite(past))) {
        r = track - past
      }
    }
    scores = if (symmetric) list(both = r) else list(lower = -r, upper = r)
    scores = lapply(scores, function(s) sort(s[is.finite(s)]))

    # Step size and integral scale from the validation scores
    ideal = apply(abs(z_cal), 2, function(v) capped_quantile(sort(v[is.finite(v)]), alpha))
    spread = if (length(ideal) > 1) diff(range(ideal)) else 1
    eta = if (method == "AC") rate else rate * spread
    KI_h = if (is.null(KI)) max(abs(z_cal[is.finite(z_cal)])) else KI
    use_I = integrate && method != "AC"

    # Starting values from the validation sample
    level = stats::setNames(rep(target, length(sides)), sides)
    P = sapply(scores, capped_quantile, a = target)
    x = stats::setNames(rep(0, length(sides)), sides)
    n_step = ncol(z_cal)
    hist = track

    n_o = length(o)
    lower = matrix(NA_real_, nrow(Fh), n_o)
    upper = lower
    centre = lower
    err = matrix(NA_real_, length(sides), n_o, dimnames = list(sides, NULL))
    used = rep(FALSE, n_o)
    tally = c(0, 0, 0)

    for (i in seq_len(n_o)) {
      t = o[i]

      # Error rates of the bands whose targets have been observed by origin t
      for (j in which(!used & o + hh <= t)) {
        used[j] = TRUE
        if (anyNA(err[, j])) {
          next
        }
        if (method == "AC") {
          level = level + eta * (target - err[, j])
        } else {
          P = P + eta * (err[, j] - target)
        }
        x = x + (err[, j] - target)
        n_step = n_step + 1
        if (method == "AcMCP" || use_scorecast) {
          z_new = (Yh[, j] - Fh[, j]) / gamma
          hist = cbind(hist, if (symmetric) abs(z_new) else z_new)
        }
      }

      # Multipliers of the lower and upper edges
      q = if (method == "AC") {
        mapply(capped_quantile, scores, level)
      } else if (use_I) {
        P + sapply(x, pid_integrator, n = n_step, KI = KI_h)
      } else {
        P
      }
      q_lower = q[[if (symmetric) "both" else "lower"]]
      q_upper = q[[if (symmetric) "both" else "upper"]]

      # AcMCP centre
      s_hat = 0
      if (method == "AcMCP") {
        s_hat = rowMeans(hist)
        if (!is.null(basis)) {
          mu = ma_intercept(as.numeric(crossprod(basis, hist)), hh)
          if (is.finite(mu)) {
            s_hat = s_hat + (mu - sum(basis * s_hat)) * basis
          }
        }
        if (use_regression && hh > 1) {
          plug = sapply(seq_len(hh - 1), function(jj) {
            k = match(t, origins_by_h[[jj]])
            if (is.null(out[[jj]]) || is.na(k)) NA_real_ else sum(basis * out[[jj]]$centre[, k])
          })
          xi_hat = cross_horizon_score(scores_pc, hh, t, plug)
          if (is.finite(xi_hat)) {
            s_hat = s_hat + (xi_hat - sum(basis * s_hat)) * basis / 2
          }
        }
      } else if (use_scorecast) {
        s_hat = score_forecast(hist, hh, scorecast)
      }

      # Band at origin t
      shift = if (symmetric) 0 else s_hat
      widen = if (symmetric) s_hat else 0
      lower[, i] = Fh[, i] + (shift - q_lower - widen) * gamma
      upper[, i] = Fh[, i] + (shift + q_upper + widen) * gamma
      centre[, i] = shift
      keep = is.finite(Yh[, i]) & !is.na(lower[, i]) & !is.na(upper[, i])
      if (any(keep)) {
        y = Yh[keep, i]
        lo = lower[keep, i]
        up = upper[keep, i]
        err[, i] = c(both = mean(y < lo | y > up), lower = mean(y < lo), upper = mean(y > up))[sides]
        score = if (all(is.finite(lo) & is.finite(up))) sum(interval_score_cells(y, lo, up, alpha)) else Inf
        tally = tally + c(sum(y >= lo & y <= up), score, length(y))
      }
    }

    out[[hh]] = list(lower = lower, upper = upper, centre = centre, origins = o, tally = tally, eta = eta)
  }

  out
}

# Learning rate of an online band with the smallest mean interval score
online_select_rate = function(method, symmetric, E_by_h, F_by_h, Y_by_h, origins_by_h, split_origin, alpha, grid,
                              integrate = TRUE, KI = NULL, regression = TRUE, scorecast = "none") {
  grid = sort(unique(grid))
  first = lapply(origins_by_h, function(o) o < split_origin)
  usable = vapply(first, function(k) sum(k) >= 2 && sum(!k) >= 1, logical(1))

  # Reference from the first half, and the run over the second half
  reference = online_reference(lapply(seq_along(E_by_h), function(hh) if (usable[hh]) E_by_h[[hh]][, first[[hh]], drop = FALSE]),
                               lapply(seq_along(E_by_h), function(hh) origins_by_h[[hh]][first[[hh]]]))
  second = lapply(seq_along(E_by_h), function(hh) !first[[hh]] & usable[hh])
  F_run = lapply(seq_along(F_by_h), function(hh) F_by_h[[hh]][, second[[hh]], drop = FALSE])
  Y_run = lapply(seq_along(Y_by_h), function(hh) Y_by_h[[hh]][, second[[hh]], drop = FALSE])
  o_run = lapply(seq_along(origins_by_h), function(hh) origins_by_h[[hh]][second[[hh]]])

  criterion = sapply(grid, function(rate) {
    run = online_bands(method, symmetric, F_run, Y_run, o_run, reference, alpha, rate, integrate, KI, regression, scorecast)
    mis = sapply(run, function(r) if (is.null(r) || r$tally[3] == 0) NA_real_ else r$tally[2] / r$tally[3])
    if (any(is.finite(mis))) mean(mis[is.finite(mis)]) else NA_real_
  })
  if (!any(is.finite(criterion))) {
    return(list(rate = 0, criterion = criterion, grid = grid))
  }
  rounded = round(criterion, 10)
  list(rate = grid[which(rounded == min(rounded, na.rm = TRUE))[1]], criterion = criterion, grid = grid)
}

# Coverage count, interval score sum and number of points for one band at one origin
band_tally = function(y, band, alpha) {
  keep = is.finite(y) & is.finite(band$lower) & is.finite(band$upper)
  if (!any(keep)) {
    return(c(0, 0, 0))
  }
  c(sum(y[keep] >= band$lower[keep] & y[keep] <= band$upper[keep]),
    sum(interval_score_cells(y[keep], band$lower[keep], band$upper[keep], alpha)),
    sum(keep))
}

band_names = c(SR = "scaled residual", NS = "non-conformity score", AC = "adaptive conformal",
               PID = "conformal PID control", AcMCP = "autocorrelated multistep conformal")
