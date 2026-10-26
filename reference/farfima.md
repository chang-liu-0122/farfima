# Fit a FARFIMA model

Fits a functional autoregressive fractionally integrated moving-average
(FARFIMA) model to the training sample and calibrates conformal
prediction bands on the validation sample.

## Usage

``` r
farfima(
  train,
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
  verbose = TRUE
)
```

## Arguments

- train:

  Training sample: a matrix of curves (grid points by time) or an `fts`
  object.

- valid:

  Validation sample, following `train`. It is needed for the conformal
  prediction bands and for `order = "cv"`.

- d:

  Memory parameter estimator for
  [`estimate_memory()`](https://chang-liu-0122.github.io/farfima/reference/estimate_memory.md)
  (`"lw"`, `"dfa"`, `"lpwn1"`, `"lpwn2"`, `"bc_lpwn1"` or `"bc_lpwn2"`),
  or a known value.

- order:

  Order selection method for
  [`select_order()`](https://chang-liu-0122.github.io/farfima/reference/select_order.md)
  (`"ic_full"`, `"ic_auto"`, `"bj"` or `"cv"`), or known orders
  `c(p, q)`. With `"cv"`, the forecasts are compared with the validation
  curves on their original scale.

- p, q:

  Candidate AR and MA orders searched by the selection method. Not used
  when `order` gives the orders. (0, 0) is excluded.

- method:

  `"farma"` for the functional Hannan-Rissanen estimator, or `"far"` for
  the FAR(1) benchmark fitted by
  [`far::far()`](https://rdrr.io/pkg/far/man/far.html). The benchmark
  ignores `order`, and with `d = 1` it is the FARFIMA(1, 1, 0)
  benchmark.

- h:

  Maximum forecast horizon of the conformal prediction bands and of
  `order = "cv"`. `valid` needs at least `h + 1` curves.

- level:

  Nominal coverage levels of the conformal prediction bands.

- band:

  Conformal prediction bands to calibrate: `"SR"` (scaled residual),
  `"NS"` (non-conformity score), `"AC"` (adaptive conformal), `"PID"`
  (conformal PID control, proportional-integral-derivative) and
  `"AcMCP"` (autocorrelated multistep conformal). See
  [`predict()`](https://chang-liu-0122.github.io/farfima/reference/predict.farfima.md).

- window:

  `"expanding"` or `"fixed"` window at each forecast origin. The fixed
  window has the length of the training sample.

- valid_eval:

  Optional curves to check the validation forecasts against, in place of
  `valid`. The forecasts are still made from `valid`. Use it when the
  model is fitted to smoothed curves but the bands should be calibrated
  on the raw curves. It must have the same size as `valid`.

- Npc:

  Number of principal components in the operator estimator.

- trim:

  Number of initial fractionally differenced curves dropped before order
  selection and estimation. Defaults to \\\min(100, \lfloor T/3
  \rfloor)\\, where \\T\\ is the number of training curves.

- far_pve:

  Proportion of variance explained by the components of the FAR(1)
  benchmark.

- memory_args:

  List of other arguments to
  [`estimate_memory()`](https://chang-liu-0122.github.io/farfima/reference/estimate_memory.md),
  such as `m_const` or `B`.

- order_args:

  List of other arguments to
  [`select_order()`](https://chang-liu-0122.github.io/farfima/reference/select_order.md),
  such as `n_lags` or `alpha`.

- eta_grid:

  Candidate learning rates of the AC, PID and AcMCP conformal prediction
  bands. One rate is chosen for each band and level, by the mean
  interval score on the second half of the validation sample.

- zeta_grid:

  Candidate multipliers of the SR band.

- symmetric:

  If `TRUE`, the AC and PID bands are symmetric about the point
  forecast. If `FALSE`, their lower and upper edges are calibrated
  separately, each at level \\\alpha/2\\. The SR and NS bands are always
  symmetric and the AcMCP band is always asymmetric.

- integrate:

  Whether the PID and AcMCP bands include the integral term, which
  widens the band after repeated under-coverage and narrows it after
  repeated over-coverage. In the extreme, the band becomes infinite or
  empty.

- KI:

  Scale of the integral term. Defaults to the largest absolute
  standardised validation residual at each horizon.

- scorecast:

  Scorecaster of the PID band: `"none"`, `"naive"` (the latest observed
  score curve) or `"theta"` (the theta method: exponential smoothing of
  the past scores plus half their linear trend).

- acmcp_regression:

  Whether the centre of the AcMCP band also uses the regression on the
  scores at shorter horizons. If `FALSE`, only the MA(h - 1) forecast is
  used.

- seed:

  Optional seed for the BC-LPWN bootstrap and the white-noise test of
  `order = "bj"`.

- verbose:

  Whether to show progress messages and a progress bar.

## Value

An object of class `farfima`, a list with

- `d` and `order`: the memory parameter and `c(p = , q = )`;

- `phi`, `theta`: the AR and MA operators (grid by grid by lag);

- `fitted`, `fitted_differenced`, `residuals`: the in-sample fit on the
  scale of the data and on the fractionally differenced scale, and the
  residuals, also given by
  [`fitted()`](https://chang-liu-0122.github.io/farfima/reference/fitted.farfima.md)
  and
  [`residuals()`](https://chang-liu-0122.github.io/farfima/reference/fitted.farfima.md)
  (`NULL` for the FAR(1) benchmark);

- `memory`, `order_selection`: the
  [`estimate_memory()`](https://chang-liu-0122.github.io/farfima/reference/estimate_memory.md)
  and
  [`select_order()`](https://chang-liu-0122.github.io/farfima/reference/select_order.md)
  results, or `NULL` when `d` or the orders are given;

- `operators`: the
  [`estimate_operators()`](https://chang-liu-0122.github.io/farfima/reference/estimate_operators.md)
  result, or `NULL` for the FAR(1) benchmark;

- `calibration`: a table of the calibrated conformal prediction bands by
  `level` and `h`: `SR_zeta` (SR multiplier), `NS_quantile` (NS
  quantile), `AC_eta`, `PID_eta` and `AcMCP_eta` (learning rates, for
  the bands chosen) and `n_origin` (number of validation origins).
  `NULL` without a validation sample.

[`print()`](https://rdrr.io/r/base/print.html) and
[`summary()`](https://rdrr.io/r/base/summary.html) show the fit,
[`coef()`](https://rdrr.io/r/stats/coef.html) gives `d`, `phi` and
`theta`, and
[`plot()`](https://chang-liu-0122.github.io/farfima/reference/plot.farfima.md)
draws the fitted or residual curves.

## Details

The fit has three steps. The memory parameter is estimated by
[`estimate_memory()`](https://chang-liu-0122.github.io/farfima/reference/estimate_memory.md).
The training curves are then centred and fractionally differenced, and
the first `trim` differenced curves are dropped, which gives
short-memory curves. Their FARMA orders are selected as in
[`select_order()`](https://chang-liu-0122.github.io/farfima/reference/select_order.md),
and their operators are estimated by the functional Hannan-Rissanen
method, as in
[`estimate_operators()`](https://chang-liu-0122.github.io/farfima/reference/estimate_operators.md).

With a validation sample, forecasts up to `h` steps ahead are made from
every validation origin, and their errors calibrate the conformal
prediction bands for each horizon and level. A validation origin whose
fit fails is skipped with a warning. A warning is also given when the
fitted operators are unstable.

## References

Kuenzer, T. (2024). Estimation of functional ARMA models. *Bernoulli*.

Li, D., Robinson, P. M. and Shang, H. L. (2020). Long-range dependent
curve time series. *Journal of the American Statistical Association*.

## See also

[`predict()`](https://chang-liu-0122.github.io/farfima/reference/predict.farfima.md)
to forecast from the fit;
[`estimate_memory()`](https://chang-liu-0122.github.io/farfima/reference/estimate_memory.md),
[`select_order()`](https://chang-liu-0122.github.io/farfima/reference/select_order.md)
and
[`estimate_operators()`](https://chang-liu-0122.github.io/farfima/reference/estimate_operators.md)
for the steps on their own.

## Examples

``` r
sim = sim_farfima(n = 200, d = 0.3, p = 1, q = 0, no_grid = 51, seed = 1)
blocks = split_curves(sim)
fit = farfima(blocks$train, blocks$valid, d = "lpwn2", order = "ic_full",
               p = 0:2, q = 0:2, h = 5, trim = 30)
#> Selecting (p, q) by functional IC over 8 candidates
#> Calibrating the SR, NS, AC bands on 40 validation origins (36 at h = 5)
fit
#> FARFIMA(p = 1, d = 0.2924, q = 0)
#>   Data:      51 grid points x 120 training curves, 40 validation curves
#>   Memory:    LPWN2 on the leading score process of the empirical covariance operator (Covar), bandwidth exponent 0.75
#>   Orders:    functional IC (exhaustive search), 8 candidates
#>   Operators: functional Hannan-Rissanen estimator with 1 principal component, first 30 fractionally differenced curves trimmed
#>   Bands:     SR, NS, AC at 80%, 95% for h = 1..5, calibrated on 40 validation curves (expanding window)
```
