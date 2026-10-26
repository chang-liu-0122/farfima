# Split curves into training, validation and test samples

Splits the curves in time order.

## Usage

``` r
split_curves(x, train = 0.6, valid = 0.2)
```

## Arguments

- x:

  A matrix of curves (grid points by time) or an `fts` object.

- train, valid:

  Proportions of the curves in the training and validation samples,
  rounded down. The remaining curves form the test sample.

## Value

A list with `train`, `valid` and `test`. `test` is `NULL` when
`train + valid = 1`.

## See also

[`farfima()`](https://chang-liu-0122.github.io/farfima/reference/farfima.md)

## Examples

``` r
x = sim_farfima(n = 100, d = 0.2, p = 1, q = 0, no_grid = 21, seed = 1)
blocks = split_curves(x, train = 0.6, valid = 0.2)
sapply(blocks, ncol)
#> train valid  test 
#>    60    20    20 
```
