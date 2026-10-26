# Estimate the memory parameter

Estimates the memory parameter \\d\\ from the leading score process, the
centred curves projected onto the leading eigenfunction of a
covariance-type operator. If the estimate is near the upper end of its
search interval, the score process is integer differenced and \\d\\ is
estimated again, until the estimate falls below the upper end. The
number of differences is then added back to the estimate.

## Usage

``` r
estimate_memory(
  x,
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
  verbose = FALSE
)
```

## Arguments

- x:

  A matrix of curves (grid points by time).

- method:

  Estimator: `"lw"` (local Whittle), `"dfa"` (detrended fluctuation
  analysis), `"lpwn1"` or `"lpwn2"` (local polynomial Whittle with
  noise, with one or two polynomial correction terms), or `"bc_lpwn1"`
  or `"bc_lpwn2"` (BC-LPWN1 or BC-LPWN2, the LPWN estimate
  bias-corrected by the prefiltered sieve bootstrap).

- cov_type:

  Operator of the FPCA that gives the leading eigenfunction: `"sample"`
  (empirical covariance operator), `"lrc_RS"` (long-run covariance
  operator with an adaptive bandwidth) or `"lrc_LRS"` (long-run
  covariance operator in sandwich form).

- m_const:

  Bandwidth exponent \\a\\, so that \\\lfloor 1 + T^a \rfloor\\
  frequencies are used.

- B, burn_in:

  Number of bootstrap samples and burn-in length (BC-LPWN).

- R_noise:

  Number of polynomial terms for the noise in LPWN.

- tuning_const:

  Constant \\c\\. The score process is differenced while the estimate is
  at least \\1 - c / \log m\\.

- cov_band_const:

  Bandwidth constant for `cov_type = "lrc_LRS"`.

- dfa_levels, dfa_minnpts, dfa_cutoff:

  Number of block sizes, minimum number of points per block and range of
  block sizes in the DFA regression.

- seed:

  Optional seed for the bootstrap.

- verbose:

  Whether to show progress.

## Value

An object of class `farfima_memory` with the estimate `d` and the
leading score process `scores`. BC-LPWN also returns the uncorrected
estimate `d_initial`, the `bias` and the `bootstrap` estimates.

## Details

The local Whittle and LPWN code is adapted from the LongMemoryTS package
(Leschinski, Voges and Wenger, 2019), and the DFA code from the fArma
package (Wuertz, Setz and Chalabi, 2017).

## References

Robinson, P. M. (1995). Gaussian semiparametric estimation of long range
dependence. *The Annals of Statistics*.

Peng, C.-K., Buldyrev, S. V., Havlin, S., Simons, M., Stanley, H. E. and
Goldberger, A. L. (1994). Mosaic organization of DNA nucleotides.
*Physical Review E*.

Frederiksen, P., Nielsen, F. S. and Nielsen, M. O. (2012). Local
polynomial Whittle estimation of perturbed fractional processes.
*Journal of Econometrics*.

Paparoditis, E. (2018). Sieve bootstrap for functional time series. *The
Annals of Statistics*.

Li, D., Robinson, P. M. and Shang, H. L. (2020). Long-range dependent
curve time series. *Journal of the American Statistical Association*.

Rice, G. and Shang, H. L. (2017). A plug-in bandwidth selection
procedure for long-run covariance estimation with stationary functional
time series. *Journal of Time Series Analysis*.

Liu, C. and Shang, H. L. (2026). Bias correction of long-memory
estimator of functional time series via the prefiltered sieve bootstrap.
arXiv:2607.28177.

Leschinski, C., Voges, M. and Wenger, K. (2019). LongMemoryTS: Long
Memory Time Series. R package version 0.1.0.

Wuertz, D., Setz, T. and Chalabi, Y. (2017). fArma: Rmetrics - Modelling
ARMA Time Series Processes. R package version 3042.81.

## See also

[`frac_diff()`](https://chang-liu-0122.github.io/farfima/reference/frac_diff.md)
to difference the curves with the estimate, and
[`farfima()`](https://chang-liu-0122.github.io/farfima/reference/farfima.md)
for the whole fit.

## Examples

``` r
x = sim_farfima(n = 300, d = 0.3, p = 1, q = 0, seed = 1)
estimate_memory(x, method = "lw")
#> Memory estimate d = 0.6358
#>   LW on the leading score process of the empirical covariance operator (Covar), bandwidth exponent 0.75
estimate_memory(x, method = "lpwn2")
#> Memory estimate d = 0.5871
#>   LPWN2 on the leading score process of the empirical covariance operator (Covar), bandwidth exponent 0.75
```
