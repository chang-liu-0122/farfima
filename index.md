  

![Macquarie University](reference/figures/macquarie-logo.png)    ![The
University of Sydney](reference/figures/usyd-logo.png)

### farfima: Order Selection and Forecasting in Long-Memory Functional Time Series

Chang Liu, Han Lin Shang, Garth Tarr and Samuel Müller  
[**Package website »**](https://chang-liu-0122.github.io/farfima/)

------------------------------------------------------------------------

## Overview

farfima provides tools for modelling, forecasting and visualising
long-memory functional time series with the functional autoregressive
fractionally integrated moving-average (FARFIMA) model. The model is
fitted in three steps, each of which can also be used on its own. First,
the memory parameter is estimated from the leading principal component
scores, for example by local Whittle or by bias-corrected local
polynomial Whittle with noise. Second, the autoregressive and
moving-average orders of the fractionally differenced curves are
selected by a functional information criterion, functional Box-Jenkins
or cross-validation. Third, the operators are estimated by the
functional Hannan-Rissanen method. The point forecasts come with
conformal prediction bands that quantify their uncertainty.
Visualisation tools include improved rainbow plots and plots of the
fitted model and its forecasts.

------------------------------------------------------------------------

## Installation

``` r

# install.packages("remotes")
remotes::install_github("chang-liu-0122/farfima", build_vignettes = TRUE)
```

farfima needs R 4.1.0 or later. Installing it also installs the packages
it uses: FCVAR, far, ftsa, longmemo, rainbow, sde, vars, withr and
wwntests. The examples below and the vignette also use the demography
package.

------------------------------------------------------------------------

## Getting Started

Curves are stored as a matrix with grid points in rows and time in
columns, such as ages by years. `fts` objects from the rainbow and ftsa
packages also work. The example uses French female mortality rates from
the demography package, for the years 1901 to 2006, with ages 0 to 94
and a group for ages 95 and above.

``` r

library(farfima)
library(demography)

# French female log death rates
fr = extract.ages(extract.years(fr.mort, 1901:2006), 0:95)
x = log10(smooth.demogdata(fr, method = "mspline", b = 65, weight = TRUE)$rate$female)
plot_curves(x)

# Training, validation and test samples
blocks = split_curves(x, train = 0.6, valid = 0.2)

# Fit the model and calibrate the conformal prediction bands
fit = farfima(blocks$train, blocks$valid, h = 10, level = c(0.8, 0.95), trim = 20, seed = 1)
fit
plot(fit, type = "residuals")

# Forecast, and evaluate on the test sample
fc = predict(fit, test = blocks$test)
fc
plot(fc, h = 1, band = "AC", xlab = "Age", ylab = "log10 death rate")
```

The test evaluation reports the RMSPE and MAPE of the point forecasts,
and the ECP, CPD, MIS and WIS of the conformal prediction bands.

### Main options of `farfima()`

| Argument | Sets | Options |
|----|----|----|
| `d` | memory parameter | `"lw"`, `"dfa"`, `"lpwn1"`, `"lpwn2"`, `"bc_lpwn1"`, `"bc_lpwn2"` (default), or a known value |
| `order` | FARMA orders | `"ic_full"` (default), `"ic_auto"`, `"bj"`, `"cv"`, or known orders `c(p, q)` |
| `p`, `q` | candidate orders | e.g. `p = 0:5, q = 0:5` |
| `method` | model | `"farma"` (default), or `"far"` for the FAR(1) benchmark |
| `h` | longest forecast horizon | e.g. `h = 10` |
| `level` | levels of the conformal prediction bands | e.g. `c(0.8, 0.95)` |
| `band` | conformal prediction bands | `"SR"`, `"NS"`, `"AC"`, `"PID"`, `"AcMCP"` |
| `window` | window at each forecast origin | `"expanding"` (default) or `"fixed"` |

### Steps of the analysis

The analysis is composed of four steps.
[`farfima()`](https://chang-liu-0122.github.io/farfima/reference/farfima.md)
runs steps 1 to 3 together and calibrates the conformal prediction
bands, and [`predict()`](https://rdrr.io/r/stats/predict.html) makes the
forecasts in step 4. Steps 1 to 3 can also be run on their own, with the
functions below.

| Step | Function | Methods |
|----|----|----|
| 1\. Memory parameter | [`estimate_memory()`](https://chang-liu-0122.github.io/farfima/reference/estimate_memory.md) | local Whittle, DFA, LPWN, bias-corrected LPWN |
| 2\. FARMA orders | [`select_order()`](https://chang-liu-0122.github.io/farfima/reference/select_order.md) | functional information criterion, functional Box-Jenkins, cross-validation |
| 3\. Operators | [`estimate_operators()`](https://chang-liu-0122.github.io/farfima/reference/estimate_operators.md) | functional Hannan-Rissanen |
| 4\. Forecasts | [`predict()`](https://rdrr.io/r/stats/predict.html) | point forecasts with conformal prediction bands (SR, NS, AC, PID, AcMCP), evaluated on a test sample |

### Plots

| Code | Shows |
|----|----|
| `plot_curves(x)` | rainbow plot of the curves, or one panel per period |
| `plot(fit)` | the fitted curves, or the residual curves with `type = "residuals"` |
| `plot(fit$order_selection)` | the criterion over the candidate orders |
| `plot(fit$operators)` | heatmaps of the estimated operators |
| `plot(fc)` | forecasts with their conformal prediction bands, by curve, horizon, grid point or recent years |

------------------------------------------------------------------------

## Folders and Files

### Package files

| File | Purpose |
|----|----|
| `DESCRIPTION` | Package name, version, authors, licence and the packages farfima uses |
| `NAMESPACE` | Functions that farfima exports and imports, generated by roxygen2 |
| `README.md` | This page, which is also the home page of the package website |
| `_pkgdown.yml` | Layout of the package website |
| `.github/workflows/pkgdown.yaml` | Instructions for GitHub to rebuild the package website after each push |
| `.Rbuildignore`, `.gitignore` | Files left out of the package build and of git |
| `farfima.Rproj` | RStudio project |

### `R/`: source code

Exported functions, one file each:

| File | Purpose |
|----|----|
| `farfima.R` | Fit the FARFIMA model and calibrate the conformal prediction bands |
| `predict.farfima.R` | Point forecasts with conformal prediction bands, and their evaluation on a test sample |
| `split_curves.R` | Split the curves into training, validation and test samples |
| `plot_curves.R` | Rainbow plots of the curves, in one panel or one panel per period |
| `estimate_memory.R` | Memory parameter estimation |
| `frac_diff.R` | Fractional differencing or integration of functional time series |
| `select_order.R` | FARMA order selection |
| `estimate_operators.R` | FARMA operator estimation |
| `sim_farfima.R` | Simulation of FARFIMA(p, d, q) curves |

Methods for fits and forecasts, including visualisation:

| File | Purpose |
|----|----|
| `fitted.farfima.R`, `residuals.farfima.R`, `coef.farfima.R` | Fitted curves, residual curves and estimates of a fit |
| `plot.farfima.R` | Visualisation of the fitted and residual curves |
| `plot.farfima_forecast.R` | Visualisation of the forecasts and conformal prediction bands |
| `plot.farfima_order.R`, `plot.farfima_operators.R` | Visualisation of the order selection criterion, and heatmaps of the operators |
| `as.data.frame.farfima_forecast.R` | Forecasts and conformal prediction bands as a data frame |
| `print.*.R`, `summary.*.R` | Printed summaries of fits, forecasts and estimates |

Internal functions:

| File | Purpose |
|----|----|
| `lw.R`, `dfa.R`, `lpwn.R` | Local Whittle, DFA and LPWN estimators |
| `sieve_bootstrap.R` | Sieve bootstrap for the bias-corrected LPWN estimators |
| `covariance.R` | Covariance-type operators and the leading principal component |
| `order_selection.R`, `ic_internal.R` | Functional information criterion, functional Box-Jenkins and cross-validation |
| `farma.R` | Functional Hannan-Rissanen estimation and the stability check |
| `forecast_paths.R`, `origin.R` | Forecasts on the original scale at each forecast origin |
| `calibrate.R`, `conformal.R` | Calibration of the SR, NS, AC, PID and AcMCP conformal prediction bands |
| `simulate_kernels.R` | Brownian motion innovations and operator kernels for [`sim_farfima()`](https://chang-liu-0122.github.io/farfima/reference/sim_farfima.md) |
| `input.R`, `labels.R`, `plot_helpers.R` | Input checks, printed labels and plotting helpers |
| `farfima-package.R` | The package help page, `?farfima-package` |

### `man/`: help pages

The help pages, generated by roxygen2 from the comments in `R/`.
`man/figures/` holds the logos on this page.

### `tests/testthat/`: tests

One test file per function, such as `test-farfima.R` and
`test-predict.R`, and `helper-curves.R` with the simulated curves shared
by the tests. Run them with `devtools::test()`.

### `vignettes/`: tutorial

`farfima.Rmd` introduces each function in turn, with an example on the
French mortality data. It is the “Get started” page of the package
website, and
[`vignette("farfima")`](https://chang-liu-0122.github.io/farfima/articles/farfima.md)
opens it in R.

------------------------------------------------------------------------

## Documentation

- The package website: <https://chang-liu-0122.github.io/farfima/>
- [`vignette("farfima")`](https://chang-liu-0122.github.io/farfima/articles/farfima.md):
  a walk-through of the package
- [`?farfima`](https://chang-liu-0122.github.io/farfima/reference/farfima.md):
  the help page of the main function, with links to the others

------------------------------------------------------------------------

## Authors

| Name | Role | Affiliation | Email |
|----|----|----|----|
| Chang Liu | Author, maintainer | Macquarie University | <chang.liu45@students.mq.edu.au> |
| Han Lin Shang | Author | Macquarie University | <hanlin.shang@mq.edu.au> |
| Garth Tarr | Author | The University of Sydney | <garth.tarr@sydney.edu.au> |
| Samuel Müller | Author | Macquarie University | <samuel.muller@mq.edu.au> |

------------------------------------------------------------------------

## License

GPL-2
