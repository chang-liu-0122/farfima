# farfima

Modelling, forecasting and visualising long-memory functional time series with
the functional autoregressive fractionally integrated moving-average (FARFIMA)
model.

## Overview

The package works in four steps. `farfima()` runs the first three together
and calibrates the conformal prediction bands, and `predict()` makes the
forecasts. Steps 1 to 3 can also be run on their own; steps 2 and 3 take
short-memory curves, such as fractionally differenced curves.

| Step | Function | Methods |
|---|---|---|
| 1. Memory parameter | `estimate_memory()` | local Whittle, DFA, LPWN, bias-corrected LPWN |
| 2. FARMA orders | `select_order()` | functional information criterion, functional Box-Jenkins, cross-validation |
| 3. Operators | `estimate_operators()` | functional Hannan-Rissanen |
| 4. Forecasts | `predict()` | point forecasts with conformal prediction bands (SR, NS, AC, PID, AcMCP), evaluated on a test sample |

## Installation

```r
# install.packages("remotes")
remotes::install_github("chang-liu-0122/farfima", build_vignettes = TRUE)
```

## Getting started

Curves are stored as a matrix with grid points in rows and time in columns,
such as ages by years. `fts` objects from the rainbow and ftsa packages also
work.

```r
library(farfima)
library(demography)

# French female log death rates, ages 0 to 100, 1901 to 2006
x = log10(smooth.demogdata(extract.years(fr.mort, 1901:2006))$rate$female[1:101, ])
plot_curves(x)

# Training, validation and test samples
blocks = split_curves(x, train = 0.6, valid = 0.2)

# Fit the model and calibrate the conformal prediction bands
fit = farfima(blocks$train, blocks$valid, h = 10, level = c(0.8, 0.95), trim = 20, seed = 1)
fit

# Forecast, and evaluate on the test sample
fc = predict(fit, test = blocks$test)
fc
plot(fc, h = 1, band = "AC", xlab = "Age", ylab = "log10 death rate")
```

The test evaluation reports the RMSPE and MAPE of the point forecasts, and the
ECP, CPD, MIS and WIS of the conformal prediction bands.

## Main options of `farfima()`

| Argument | Sets | Options |
|---|---|---|
| `d` | memory parameter | `"lw"`, `"dfa"`, `"lpwn1"`, `"lpwn2"`, `"bc_lpwn1"`, `"bc_lpwn2"` (default), or a known value |
| `order` | FARMA orders | `"ic_full"` (default), `"ic_auto"`, `"bj"`, `"cv"`, or known orders `c(p, q)` |
| `p`, `q` | candidate orders | e.g. `p = 0:5, q = 0:5` |
| `method` | model | `"farma"` (default), or `"far"` for the FAR(1) benchmark |
| `h` | longest forecast horizon | e.g. `h = 10` |
| `level` | levels of the conformal prediction bands | e.g. `c(0.8, 0.95)` |
| `band` | conformal prediction bands | `"SR"`, `"NS"`, `"AC"`, `"PID"`, `"AcMCP"` |
| `window` | window at each forecast origin | `"expanding"` (default) or `"fixed"` |

## Plots

| Code | Shows |
|---|---|
| `plot_curves(x)` | rainbow plot of the curves, or one panel per period |
| `plot(fit$order_selection)` | the criterion over the candidate orders |
| `plot(fit$operators)` | heatmaps of the estimated operators |
| `plot(fc)` | forecasts with their conformal prediction bands, by curve, horizon, grid point or recent years |

## Learn more

- `vignette("farfima")`: a full walk-through of the package.
- `?farfima`: the help page of the main function, with links to the others.

## License

GPL-2
