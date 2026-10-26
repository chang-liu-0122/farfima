# Plot the fitted or residual curves of a FARFIMA fit

Rainbow plots of the in-sample fit of a model from
[`farfima()`](https://chang-liu-0122.github.io/farfima/reference/farfima.md),
coloured by time as in
[`plot_curves()`](https://chang-liu-0122.github.io/farfima/reference/plot_curves.md).
`type = "fitted"` draws the fitted curves on the scale of the data, and
`type = "residuals"` draws the residual curves. The first `max(p, q)`
curves have no fit and are left out.

## Usage

``` r
# S3 method for class 'farfima'
plot(
  x,
  type = c("fitted", "residuals"),
  argvals = NULL,
  xlab = "Grid",
  ylab = NULL,
  main = NULL,
  ...
)
```

## Arguments

- x:

  A fit from
  [`farfima()`](https://chang-liu-0122.github.io/farfima/reference/farfima.md)
  with `method = "farma"`.

- type:

  `"fitted"` or `"residuals"`.

- argvals:

  Grid values. Defaults to the row names of the curves if numeric.

- xlab, ylab, main:

  Axis labels and title.

- ...:

  Other arguments passed to
  [`graphics::matplot()`](https://rdrr.io/r/graphics/matplot.html).

## Value

`x`, invisibly.

## See also

[`fitted()`](https://chang-liu-0122.github.io/farfima/reference/fitted.farfima.md)
and
[`residuals()`](https://chang-liu-0122.github.io/farfima/reference/fitted.farfima.md)
for the curves.

## Examples

``` r
sim = sim_farfima(n = 200, d = 0.3, p = 1, q = 0, no_grid = 51, seed = 1)
fit = farfima(sim, d = 0.3, order = c(1, 0), trim = 50, verbose = FALSE)
plot(fit)

plot(fit, type = "residuals")
```
