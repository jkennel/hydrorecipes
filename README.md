
<!-- README.md is generated from README.Rmd. Please edit that file -->

# hydrorecipes

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/hydrorecipes)](https://CRAN.R-project.org/package=hydrorecipes)
[![R-CMD-check](https://github.com/jkennel/hydrorecipes/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jkennel/hydrorecipes/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The goal of hydrorecipes is to supplement the recipes package with a few
steps that can help deal with moderately sized water level datasets.

## Installation

You can install the development version of hydrorecipes from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("jkennel/hydrorecipes")
```

## Example

## Kennel 2020

This is the method of Kennel 2020 which uses a distributed lag model and
the **earthtide** package to generate synthetic wave groups.

``` r
library(hydrorecipes)
#> Loading required package: recipes
#> Loading required package: dplyr
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
#> 
#> Attaching package: 'recipes'
#> The following object is masked from 'package:stats':
#> 
#>     step
library(earthtide)

data(transducer)

# convert to numeric because step_ns doesn't handle POSIXct
transducer$datetime <- as.numeric(transducer$datetime)

# Earth tide inputs
wave_groups <- earthtide::eterna_wavegroups
wave_groups <- na.omit(wave_groups[wave_groups$time == '1 month', ])
wave_groups <- wave_groups[wave_groups$start > 0.5,]
latitude  <- 34.0
longitude <- -118.5
elevation <- 500
cutoff    <- 1e-6
catalog   <- 'ksm04'
method    <- 'volume_strain'

# create recipe 
rec <- recipe(wl~baro+datetime, transducer) |>
  step_distributed_lag(baro, knots = log_lags(15, 86400 * 2 / 120)) |>
  step_earthtide(datetime,
                 latitude = latitude,
                 longitude = longitude,
                 elevation = elevation,
                 cutoff = cutoff,
                 method = method,
                 catalog = catalog,
                 astro_update = 1,
                 wave_groups = wave_groups) |>
  step_ns(datetime, deg_free = 10) |>
  prep()

input <- rec |> bake(new_data = NULL)
```

``` r
summary(fit <- lm(wl~., input))
#> 
#> Call:
#> lm(formula = wl ~ ., data = input)
#> 
#> Residuals:
#>        Min         1Q     Median         3Q        Max 
#> -0.0014831 -0.0002350 -0.0000195  0.0002125  0.0016613 
#> 
#> Coefficients:
#>                             Estimate Std. Error   t value Pr(>|t|)    
#> (Intercept)                5.321e+00  2.254e-03  2360.723  < 2e-16 ***
#> distributed_lag_baro_0    -1.940e-01  1.153e-02   -16.821  < 2e-16 ***
#> distributed_lag_baro_1     2.381e-02  8.788e-03     2.710  0.00674 ** 
#> distributed_lag_baro_2     7.713e-03  5.264e-03     1.465  0.14291    
#> distributed_lag_baro_4     1.378e-02  2.477e-03     5.563 2.66e-08 ***
#> distributed_lag_baro_7     9.760e-03  1.227e-03     7.957 1.81e-15 ***
#> distributed_lag_baro_12    7.815e-03  5.499e-04    14.211  < 2e-16 ***
#> distributed_lag_baro_22    6.232e-03  2.342e-04    26.614  < 2e-16 ***
#> distributed_lag_baro_37    3.434e-03  9.676e-05    35.487  < 2e-16 ***
#> distributed_lag_baro_63    1.634e-03  4.092e-05    39.929  < 2e-16 ***
#> distributed_lag_baro_106   1.942e-04  1.804e-05    10.765  < 2e-16 ***
#> distributed_lag_baro_179   6.993e-05  8.150e-06     8.580  < 2e-16 ***
#> distributed_lag_baro_302  -6.871e-05  4.448e-06   -15.448  < 2e-16 ***
#> distributed_lag_baro_509  -7.733e-02  1.515e-03   -51.049  < 2e-16 ***
#> distributed_lag_baro_856   2.007e-01  3.930e-03    51.071  < 2e-16 ***
#> distributed_lag_baro_1440 -1.235e-01  2.415e-03   -51.122  < 2e-16 ***
#> earthtide_cos_1            9.110e-06  3.371e-06     2.703  0.00688 ** 
#> earthtide_sin_1           -4.254e-06  3.403e-06    -1.250  0.21128    
#> earthtide_cos_2           -8.802e-05  3.474e-06   -25.333  < 2e-16 ***
#> earthtide_sin_2            3.556e-05  3.502e-06    10.155  < 2e-16 ***
#> earthtide_cos_3           -3.120e-05  2.939e-06   -10.616  < 2e-16 ***
#> earthtide_sin_3            3.341e-05  2.851e-06    11.717  < 2e-16 ***
#> earthtide_cos_4           -1.920e-04  7.844e-06   -24.477  < 2e-16 ***
#> earthtide_sin_4           -4.311e-05  6.841e-06    -6.302 2.97e-10 ***
#> earthtide_cos_5            4.425e-06  3.283e-06     1.348  0.17768    
#> earthtide_sin_5            1.429e-05  3.283e-06     4.353 1.34e-05 ***
#> earthtide_cos_6           -1.096e-05  3.706e-06    -2.957  0.00311 ** 
#> earthtide_sin_6            8.148e-06  3.720e-06     2.190  0.02852 *  
#> earthtide_cos_7           -3.210e-06  2.070e-06    -1.551  0.12100    
#> earthtide_sin_7            1.613e-05  2.072e-06     7.786 7.10e-15 ***
#> earthtide_cos_8           -5.697e-05  2.601e-06   -21.905  < 2e-16 ***
#> earthtide_sin_8            4.886e-05  2.601e-06    18.783  < 2e-16 ***
#> earthtide_cos_9           -3.015e-04  2.694e-06  -111.936  < 2e-16 ***
#> earthtide_sin_9            2.919e-04  2.681e-06   108.876  < 2e-16 ***
#> earthtide_cos_10          -4.803e-06  3.076e-06    -1.561  0.11845    
#> earthtide_sin_10          -5.247e-05  3.059e-06   -17.151  < 2e-16 ***
#> earthtide_cos_11          -2.371e-04  6.475e-06   -36.616  < 2e-16 ***
#> earthtide_sin_11           5.427e-05  5.822e-06     9.323  < 2e-16 ***
#> earthtide_cos_12          -2.096e-06  2.531e-06    -0.828  0.40776    
#> earthtide_sin_12           2.501e-06  2.536e-06     0.986  0.32403    
#> datetime_ns_01            -2.965e-02  2.892e-05 -1024.998  < 2e-16 ***
#> datetime_ns_02            -4.064e-02  3.793e-05 -1071.515  < 2e-16 ***
#> datetime_ns_03            -6.178e-02  3.313e-05 -1864.536  < 2e-16 ***
#> datetime_ns_04            -7.809e-02  3.445e-05 -2266.381  < 2e-16 ***
#> datetime_ns_05            -9.233e-02  3.317e-05 -2783.744  < 2e-16 ***
#> datetime_ns_06            -1.100e-01  3.553e-05 -3094.888  < 2e-16 ***
#> datetime_ns_07            -1.264e-01  3.354e-05 -3768.018  < 2e-16 ***
#> datetime_ns_08            -1.355e-01  2.239e-05 -6052.921  < 2e-16 ***
#> datetime_ns_09            -1.639e-01  7.196e-05 -2278.102  < 2e-16 ***
#> datetime_ns_10            -1.499e-01  1.563e-05 -9586.010  < 2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error: 0.0003564 on 35231 degrees of freedom
#>   (1440 observations deleted due to missingness)
#> Multiple R-squared:  0.9999, Adjusted R-squared:  0.9999 
#> F-statistic: 1.067e+07 on 49 and 35231 DF,  p-value: < 2.2e-16
```
