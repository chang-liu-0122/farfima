# Forecasts and conformal prediction bands as a data frame

One row per grid point, horizon, band and level, for saving the results
or for your own plots and summaries.

## Usage

``` r
# S3 method for class 'farfima_forecast'
as.data.frame(x, row.names = NULL, optional = FALSE, test = FALSE, ...)
```

## Arguments

- x:

  A forecast from
  [`predict()`](https://chang-liu-0122.github.io/farfima/reference/predict.farfima.md).

- row.names, optional:

  Not used; they come from
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html).

- test:

  `FALSE` for the forecasts after the last validation or test curve;
  `TRUE` for the forecasts of the test curves, with the observed values.

- ...:

  Not used.

## Value

A data frame with one row per grid point, horizon, band and level:

- `band`, `level`, `h`: the band type, its nominal level and the
  horizon;

- `target`: the time being forecast;

- `grid`: the grid point;

- `point_forecast`, `lower`, `upper`: the point forecast and the edges
  of the conformal prediction band at that grid point.

With `test = TRUE`, the rows cover the forecasts of the test curves,
with three more columns:

- `origin`: the time the forecast was made from;

- `target_test_curve`: the position of the target in the test sample;

- `observed`: the observed value at that grid point.

## See also

[`predict()`](https://chang-liu-0122.github.io/farfima/reference/predict.farfima.md)

## Examples

``` r
sim = sim_farfima(n = 200, d = 0.3, p = 1, q = 0, no_grid = 51, seed = 1)
blocks = split_curves(sim)
fit = farfima(blocks$train, blocks$valid, d = 0.3, order = c(1, 0), h = 3,
               trim = 50, verbose = FALSE)
fc = predict(fit, test = blocks$test, verbose = FALSE)
head(as.data.frame(fc))
#>   band level h target grid point_forecast     lower     upper
#> 1   SR   0.8 1     h1    1      -1.333967 -1.366452 -1.301482
#> 2   SR   0.8 1     h1    2      -1.427652 -1.624740 -1.230565
#> 3   SR   0.8 1     h1    3      -1.456614 -1.740898 -1.172330
#> 4   SR   0.8 1     h1    4      -1.497945 -1.863685 -1.132206
#> 5   SR   0.8 1     h1    5      -1.444699 -1.867464 -1.021935
#> 6   SR   0.8 1     h1    6      -1.448035 -1.887083 -1.008986
d = as.data.frame(fc, test = TRUE)
d = subset(d, band == "AC" & level == 0.8)
tapply(d$observed >= d$lower & d$observed <= d$upper, d$h, mean)
#>         1         2         3 
#> 0.7906863 0.7702363 0.7389061 
```
