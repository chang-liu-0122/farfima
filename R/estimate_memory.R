#' Estimate the memory parameter
#'
#' Estimates the memory parameter \eqn{d} from the leading score process, the
#' centred curves projected onto the leading eigenfunction of a
#' covariance-type operator. If the estimate is near the upper end of its
#' search interval, the score process is integer differenced and \eqn{d} is
#' estimated again, until the estimate falls below the upper end. The number
#' of differences is then added back to the estimate.
#'
#' The local Whittle and LPWN code is adapted from the LongMemoryTS package
#' (Leschinski, Voges and Wenger, 2019), and the DFA code from the fArma
#' package (Wuertz, Setz and Chalabi, 2017).
#'
#' @param x A matrix of curves (grid points by time).
#' @param method Estimator: `"lw"` (local Whittle), `"dfa"` (detrended
#'   fluctuation analysis), `"lpwn1"` or `"lpwn2"` (local polynomial Whittle
#'   with noise, with one or two polynomial correction terms), or `"bc_lpwn1"`
#'   or `"bc_lpwn2"` (BC-LPWN1 or BC-LPWN2, the LPWN estimate bias-corrected by
#'   the prefiltered sieve bootstrap).
#' @param cov_type Operator of the FPCA that gives the leading eigenfunction:
#'   `"sample"` (empirical covariance operator), `"lrc_RS"` (long-run
#'   covariance operator with an adaptive bandwidth) or `"lrc_LRS"` (long-run
#'   covariance operator in sandwich form).
#' @param m_const Bandwidth exponent \eqn{a}, so that \eqn{\lfloor 1 + T^a
#'   \rfloor} frequencies are used.
#' @param B,burn_in Number of bootstrap samples and burn-in length (BC-LPWN).
#' @param R_noise Number of polynomial terms for the noise in LPWN.
#' @param tuning_const Constant \eqn{c}. The score process is differenced while
#'   the estimate is at least \eqn{1 - c / \log m}.
#' @param cov_band_const Bandwidth constant for `cov_type = "lrc_LRS"`.
#' @param dfa_levels,dfa_minnpts,dfa_cutoff Number of block sizes, minimum
#'   number of points per block and range of block sizes in the DFA
#'   regression.
#' @param seed Optional seed for the bootstrap.
#' @param verbose Whether to show progress.
#' @return An object of class `farfima_memory` with the estimate `d` and the
#'   leading score process `scores`. BC-LPWN also returns the uncorrected
#'   estimate `d_initial`, the `bias` and the `bootstrap` estimates.
#' @seealso [frac_diff()] to difference the curves with the estimate, and
#'   [farfima()] for the whole fit.
#' @references
#' Robinson, P. M. (1995). Gaussian semiparametric estimation of long range
#' dependence. *The Annals of Statistics*.
#'
#' Peng, C.-K., Buldyrev, S. V., Havlin, S., Simons, M., Stanley, H. E. and
#' Goldberger, A. L. (1994). Mosaic organization of DNA nucleotides.
#' *Physical Review E*.
#'
#' Frederiksen, P., Nielsen, F. S. and Nielsen, M. O. (2012). Local polynomial
#' Whittle estimation of perturbed fractional processes. *Journal of
#' Econometrics*.
#'
#' Paparoditis, E. (2018). Sieve bootstrap for functional time series. *The
#' Annals of Statistics*.
#'
#' Li, D., Robinson, P. M. and Shang, H. L. (2020). Long-range dependent curve
#' time series. *Journal of the American Statistical Association*.
#'
#' Rice, G. and Shang, H. L. (2017). A plug-in bandwidth selection procedure
#' for long-run covariance estimation with stationary functional time series.
#' *Journal of Time Series Analysis*.
#'
#' Liu, C. and Shang, H. L. (2026). Bias correction of long-memory estimator
#' of functional time series via the prefiltered sieve bootstrap.
#' arXiv:2607.28177.
#'
#' Leschinski, C., Voges, M. and Wenger, K. (2019). LongMemoryTS: Long Memory
#' Time Series. R package version 0.1.0.
#'
#' Wuertz, D., Setz, T. and Chalabi, Y. (2017). fArma: Rmetrics - Modelling
#' ARMA Time Series Processes. R package version 3042.81.
#' @examples
#' x = sim_farfima(n = 300, d = 0.3, p = 1, q = 0, seed = 1)
#' estimate_memory(x, method = "lw")
#' estimate_memory(x, method = "lpwn2")
#' @export
estimate_memory = function(x,
                           method = c("bc_lpwn2", "bc_lpwn1", "lpwn2", "lpwn1", "lw", "dfa"),
                           cov_type = c("sample", "lrc_RS", "lrc_LRS"),
                           m_const = 0.75,
                           B = 100,
                           burn_in = 200,
                           R_noise = 0,
                           tuning_const = 1,
                           cov_band_const = 1,
                           dfa_levels = 50,
                           dfa_minnpts = 3,
                           dfa_cutoff = 10^c(0.7, 2.5),
                           seed = NULL,
                           verbose = FALSE) {
  method = match.arg(method)
  cov_type = match.arg(cov_type)
  x = as_curves(x, "x")
  if (!is.null(seed)) {
    withr::local_seed(seed)
  }

  # Leading score process of the centred curves
  xc = sweep(x, 1, rowMeans(x), "-")
  lead = leading_pc(xc, get_covariance(xc, type = cov_type, band_const = cov_band_const))
  R_short = if (method %in% c("lpwn1", "bc_lpwn1")) 1 else 2

  d_initial = as.numeric(switch(
    method,
    lw = local_whittle_general(lead$scores, m_const = m_const, tuning_const = tuning_const),
    dfa = peng_general(lead$scores, levels = dfa_levels, minnpts = dfa_minnpts,
                       cut.off = dfa_cutoff, tuning_const = tuning_const, m_const = m_const),
    LPWN_general(lead$scores, m_const = m_const, R_short = R_short, R_noise = R_noise,
                 tuning_const = tuning_const)
  ))

  # LPWN estimate at the boundary of its search interval?
  at_edge = function(v) {
    r = v %% 1
    is.finite(v) & (abs(r - 0.01) < 1e-6 | abs(r - 0.99) < 1e-6)
  }
  is_lpwn = method %in% c("lpwn1", "lpwn2", "bc_lpwn1", "bc_lpwn2")
  is_bc = method %in% c("bc_lpwn1", "bc_lpwn2")
  at_bound = is_lpwn && at_edge(d_initial)
  boot_at_bound = NA_real_
  d = d_initial
  bias = NA_real_
  boot_sd = NA_real_
  boot = NULL

  # BC-LPWN: correct the bias with LPWN estimates from sieve bootstrap samples
  if (is_bc) {
    x_diff = t(FCVAR::FracDiff(t(xc), d = d_initial))
    x_boot = sieve_bootstrap(fun_dat = t(x_diff), CPV_percent = 0.95, VAR_type = "none",
                             B = B, burn_in = burn_in)$X_boot
    boot = rep(NA_real_, B)
    if (verbose) {
      message(sprintf("Estimating d by BC-LPWN%d with %d bootstrap samples", R_short, B))
    }
    bar = progress_bar(verbose, B)

    for (b in seq_len(B)) {
      xb = t(FCVAR::FracDiff(t(x_boot[, , b]), d = -d_initial))
      lead_b = leading_pc(xb, get_covariance(xb, type = cov_type, band_const = cov_band_const))
      boot[b] = tryCatch(
        as.numeric(LPWN_general(lead_b$scores, m_const = m_const, R_short = R_short,
                                R_noise = R_noise, tuning_const = tuning_const)),
        error = function(e) NA_real_
      )
      if (!is.null(bar)) utils::setTxtProgressBar(bar, b)
    }
    if (!is.null(bar)) close(bar)

    boot_at_bound = mean(at_edge(boot[is.finite(boot)]))
    bias = mean(boot, na.rm = TRUE) - d_initial
    boot_sd = sd(boot - d_initial, na.rm = TRUE)
    d = d_initial - bias
  }

  admissible = is.finite(d) && (!is_bc || d >= 0)

  out = list(d = d, scores = lead$scores)
  if (is_bc) {
    out$d_initial = d_initial
    out$bias = bias
    out$bootstrap = boot
  }
  structure(out, class = "farfima_memory", method = method, cov_type = cov_type, m_const = m_const,
            sd = boot_sd, eigenfunction = lead$psi1,
            diagnostics = list(admissible = admissible, at_bound = at_bound, boot_at_bound = boot_at_bound))
}
