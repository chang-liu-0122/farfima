# Fractional differencing

Applies \\(1 - B)^d\\ at every grid point, where \\B\\ is the backshift
operator. A negative `d` integrates, so `frac_diff(frac_diff(x, d), -d)`
returns `x`.

## Usage

``` r
frac_diff(x, d)
```

## Arguments

- x:

  A matrix of curves (grid points by time).

- d:

  Memory parameter, the order of differencing. A negative `d`
  fractionally integrates the curves.

## Value

The fractionally differenced curves, or the fractionally integrated
curves when `d` is negative, as a matrix of the same size as `x`.

## See also

[`estimate_memory()`](https://chang-liu-0122.github.io/farfima/reference/estimate_memory.md)

## Examples

``` r
x = sim_farfima(n = 100, d = 0.3, p = 1, q = 0, seed = 1)
y = frac_diff(x, 0.3)
max(abs(frac_diff(y, -0.3) - x))
#> [1] 4.884981e-15
```
