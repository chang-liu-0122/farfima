# Fitted curves and residuals of a FARFIMA fit

The in-sample fit of a model from
[`farfima()`](https://chang-liu-0122.github.io/farfima/reference/farfima.md).
Each fitted curve is the one-step-ahead forecast of a training curve
from the curves before it, and each residual curve is the training curve
minus its fitted curve. They cover the training curves left after
trimming; the first `max(p, q)` columns are `NA`, because the model
needs that many earlier curves.

## Usage

``` r
# S3 method for class 'farfima'
fitted(object, scale = c("original", "differenced"), ...)

# S3 method for class 'farfima'
residuals(object, ...)
```

## Arguments

- object:

  A fit from
  [`farfima()`](https://chang-liu-0122.github.io/farfima/reference/farfima.md)
  with `method = "farma"`.

- scale:

  `"original"` or `"differenced"`.

- ...:

  Not used.

## Value

[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) returns the
fitted curves and
[`residuals()`](https://rdrr.io/r/stats/residuals.html) the residual
curves, each as a matrix with grid points in rows and training curves in
columns.

## Details

`fitted(fit)` is on the scale of the data, so
`fitted(fit) + residuals(fit)` gives the training curves.
`fitted(fit, scale = "differenced")` is the FARMA fit of the
fractionally differenced curves. The residuals are the same on both
scales.

## See also

[`farfima()`](https://chang-liu-0122.github.io/farfima/reference/farfima.md);
[`plot()`](https://chang-liu-0122.github.io/farfima/reference/plot.farfima.md)
to draw the curves.

## Examples

``` r
sim = sim_farfima(n = 200, d = 0.3, p = 1, q = 0, seed = 1)
fit = farfima(sim, d = 0.3, order = c(1, 0), trim = 50, verbose = FALSE)
dim(fitted(fit))
#> [1] 101 150
range(fitted(fit) + residuals(fit) - sim[, 51:200], na.rm = TRUE)
#> [1] -4.440892e-16  4.440892e-16
```
