# Name of the memory estimator and covariance-type operator for printing
memory_label = function(m) {
  estimators = c(bc_lpwn2 = "BC-LPWN2", bc_lpwn1 = "BC-LPWN1", lpwn2 = "LPWN2", lpwn1 = "LPWN1",
                 lw = "LW", dfa = "DFA")
  operators = c(sample = "empirical covariance operator (Covar)",
                lrc_RS = "long-run covariance operator (RS)",
                lrc_LRS = "long-run covariance operator in sandwich form (LRS)")
  paste0(estimators[[attr(m, "method")]], " on the leading score process of the ",
         operators[[attr(m, "cov_type")]], ", bandwidth exponent ", attr(m, "m_const"))
}

# Name of the order selection method for printing
order_label = function(o) {
  methods = c(ic_full = "functional IC (exhaustive search)", ic_auto = "auto functional IC",
              bj = "functional Box-Jenkins", cv = "forecast-based CV")
  method = attr(o, "method")
  extra = switch(method,
                 bj = paste0(", ", attr(o, "n_lags"), " lags"),
                 cv = paste0(", h = 1..", attr(o, "h"), ", ", attr(o, "window"), " window"),
                 "")
  paste0(methods[[method]], extra, ", ", attr(o, "n_candidates"), " candidates")
}

# Notes on LPWN estimates at the boundary of the search interval
memory_notes = function(m, indent) {
  notes = character(0)
  dg = attr(m, "diagnostics")
  if (isTRUE(dg$at_bound)) {
    notes = c(notes, paste0(indent, "note: the LPWN estimate is at the boundary of its search interval\n"))
  }
  if (isTRUE(dg$boot_at_bound > 0.25)) {
    notes = c(notes, paste0(indent, "note: ", round(100 * dg$boot_at_bound), "% of the bootstrap estimates are at the boundary\n"))
  }
  notes
}

# Nominal level as a label such as "80%"
level_label = function(level) {
  paste0(format(100 * level, trim = TRUE), "%")
}

# Time labels of the forecast curves
future_labels = function(time, h) {
  t_num = suppressWarnings(as.numeric(time))
  steps = diff(t_num)
  if (length(t_num) > 1 && all(is.finite(t_num)) && all(abs(steps - steps[1]) < 1e-8)) {
    format(t_num[length(t_num)] + steps[1] * seq_len(h), trim = TRUE)
  } else {
    paste0("h", seq_len(h))
  }
}
