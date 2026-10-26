# Simulate a FARFIMA(p, d, q) process

Simulates `n` curves from a FARFIMA(p, d, q) model on `no_grid` equally
spaced points of \[0, 1\]. The innovations are Brownian motion curves.
They are fractionally integrated with memory parameter `d`, passed
through the MA operators and then through the AR recursion. The first
100 curves are dropped as burn-in.

## Usage

``` r
sim_farfima(
  n = 500,
  d = 0.1,
  p = 1,
  q = 0,
  no_grid = 101,
  seed = NULL,
  ar_strength = 0.6,
  ma_strength = 0.6,
  lag_strength = c("equal", "decay")
)
```

## Arguments

- n:

  Number of curves.

- d:

  Memory parameter.

- p, q:

  AR and MA orders, not both zero.

- no_grid:

  Number of grid points.

- seed:

  Optional seed.

- ar_strength, ma_strength:

  Operator strengths \\\psi\\ and \\\nu\\.

- lag_strength:

  `"equal"` splits the strength equally across lags, and `"decay"` lets
  it decay geometrically.

## Value

A matrix of curves, `no_grid` by `n`. The innovations are kept in
`attr(, "innovations")`.

## Details

The AR kernel is \\\exp\\(u^2 + v^2)/2\\\\ and the MA kernel is
\\\min(u, v)\\. Each is scaled to operator norm one and multiplied by
`ar_strength` or `ma_strength`.

## See also

[`split_curves()`](https://chang-liu-0122.github.io/farfima/reference/split_curves.md),
[`farfima()`](https://chang-liu-0122.github.io/farfima/reference/farfima.md)

## Examples

``` r
sim = sim_farfima(n = 200, d = 0.2, p = 1, q = 1, seed = 1)
dim(sim)
#> [1] 101 200
```
