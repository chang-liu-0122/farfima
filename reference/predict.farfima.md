# Forecast from a FARFIMA fit

Point forecasts and conformal prediction bands up to `h` steps beyond
the data. With a test sample, forecasts are also made from every test
origin and evaluated.

## Usage

``` r
# S3 method for class 'farfima'
predict(
  object,
  test = NULL,
  h = attr(object, "settings")$h,
  band = attr(object, "settings")$band,
  test_eval = NULL,
  verbose = TRUE,
  ...
)
```

## Arguments

- object:

  A fit from
  [`farfima()`](https://chang-liu-0122.github.io/farfima/reference/farfima.md).

- test:

  Optional test sample, following the validation sample.

- h:

  Maximum forecast horizon, at most the `h` of the fit.

- band:

  Conformal prediction bands to return, from those calibrated in
  [`farfima()`](https://chang-liu-0122.github.io/farfima/reference/farfima.md).

- test_eval:

  Optional curves to check the test forecasts against, in place of
  `test`. The forecasts are still made from `test`.

- verbose:

  Whether to show progress messages and a progress bar.

- ...:

  Not used.

## Value

An object of class `farfima_forecast`, a list with

- `point_forecast`: the forecast curves beyond the data (grid by `h`);

- `bands`: the conformal prediction bands, with the `lower` and `upper`
  matrices (grid by `h`) for each band and level (such as `"80%"`);

- `evaluation`, with a test sample: the tables `point` (RMSPE, the root
  mean squared prediction error, and MAPE, the mean absolute prediction
  error), `interval` (ECP, CPD and MIS of the bands) and `wis` (weighted
  interval score).

## Details

The conformal prediction bands are those chosen in
[`farfima()`](https://chang-liu-0122.github.io/farfima/reference/farfima.md):

- `"SR"`: scaled residual band, the forecast plus and minus a multiplier
  of the point-wise scale function;

- `"NS"`: non-conformity score band, with the width set by a quantile of
  the standardised non-conformity scores;

- `"AC"`: adaptive conformal band, an NS band whose quantile level is
  updated over the test sample;

- `"PID"`: conformal PID control band, whose multiplier is updated by a
  proportional term, an integral term and an optional scorecaster;

- `"AcMCP"`: autocorrelated multistep conformal band, an asymmetric PID
  band shifted by a forecast of the standardised residual curve.

The AC, PID and AcMCP bands are updated from the error rates of earlier
bands, once their targets are observed.

In the test evaluation, the memory parameter and the orders stay fixed,
and the operators and the mean are re-estimated at each origin. An
origin whose fit fails is skipped with a warning. The forecast beyond
the data starts after the last test curve, or after the last validation
curve when there is no test sample.

## References

Gibbs, I. and Candes, E. J. (2021). Adaptive conformal inference under
distribution shift. *Advances in Neural Information Processing Systems*.

Angelopoulos, A. N., Candes, E. J. and Tibshirani, R. J. (2023).
Conformal PID control for time series prediction. *Advances in Neural
Information Processing Systems*.

Wang, X. and Hyndman, R. J. (2024). Online conformal inference for
multi-step time series forecasting. arXiv:2410.13115.

Wang, X. and Girolimetto, D. (2026). conformalForecast: an R package for
sequential conformal prediction in time series. arXiv:2610.21648.

Bracher, J., Ray, E. L., Gneiting, T. and Reich, N. G. (2021).
Evaluating epidemic forecasts in an interval format. *PLOS Computational
Biology*.

## See also

[`plot()`](https://chang-liu-0122.github.io/farfima/reference/plot.farfima_forecast.md)
and
[`as.data.frame()`](https://chang-liu-0122.github.io/farfima/reference/as.data.frame.farfima_forecast.md)
for the results.

## Examples

``` r
sim = sim_farfima(n = 200, d = 0.3, p = 1, q = 0, no_grid = 51, seed = 1)
blocks = split_curves(sim)
fit = farfima(blocks$train, blocks$valid, d = 0.3, order = c(1, 0), h = 5, trim = 30)
#> Calibrating the SR, NS, AC bands on 40 validation origins (36 at h = 5)
fc = predict(fit, test = blocks$test)
#> Evaluating on 40 test origins
fc
#> Forecasts from FARFIMA(p = 1, d = 0.3, q = 0)
#>   Horizons:  h = 1..5 after curve 200 (expanding window)
#>   Bands:     SR (scaled residual), NS (non-conformity score), AC (adaptive conformal)
#>   Levels:    80%, 95%
#> 
#> Test evaluation over 40 test curves, averaged over horizons:
#>   RMSPE = 0.8063, MAPE = 0.6193
#>  band level    ECP     CPD   MIS
#>    SR   80% 0.7798 0.04782 2.851
#>    NS   80% 0.7762 0.04989 2.853
#>    AC   80% 0.7463 0.05368 2.816
#>    SR   95% 0.9043 0.04809 4.009
#>    NS   95% 0.9040 0.04851 4.013
#>    AC   95% 0.9232 0.02679 3.753
#>   WIS: SR = 0.2780, NS = 0.2781, AC = 0.2740
```
