
<!-- README.md is generated from README.Rmd. Please edit that file -->

# hydrorecipes

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/hydrorecipes)](https://CRAN.R-project.org/package=hydrorecipes)
[![R-CMD-check](https://github.com/jkennel/hydrorecipes/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jkennel/hydrorecipes/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The goal of hydrorecipes is to supplement the [recipes
package](https://recipes.tidymodels.org) with a few steps that can help
deal with moderately sized water level datasets. These were developed
primarily with regression convolution in mind. The following steps are
currently available:

-   *step_lead_lag* is a more flexible version of *step_lag* from the
    [recipes package](https://recipes.tidymodels.org). Values can be
    negative which indicates the vector is lead. In addition, subsetting
    can be done on the dataset during this process when dealing with
    very large datasets leading to performance gains without sacrificing
    temporal accuracy.
-   *step_distributed_lag* is a distributed lag approach for modelling
    the response in a flexible yet concise manner. This is useful when
    you have a long maximum lag or large datasets.
-   *step_earthtide* uses the [earthtide
    package](https://cran.rstudio.com/web/packages/earthtide/index.html)
    to model the synthetic Earth tide given locations and times. This
    can provide a single Earth tide curve or a set of harmonics that can
    be used in regression models.

## Installation

You can install the development version of hydrorecipes from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("jkennel/hydrorecipes")
```

## Example

This is the method of [Kennel 2020](http://hdl.handle.net/10214/17890)
which uses a distributed lag model and the [earthtide
package](https://cran.rstudio.com/web/packages/earthtide/index.html) to
generate synthetic wave groups. A \~1.5 month dataset of water and
barometric pressure having a monitoring frequency of 2 minutes is
presented below. The barometric response is modeled over two days using
a distributed lag model with 12 regressor terms. The *knots* are
logarithmically separated over two days to accurately capture early and
late time responses which can be caused by different physical
mechanisms.

``` r
library(hydrorecipes)
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
  step_ns(datetime, deg_free = 12) |>
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
#> -1.596e-03 -2.392e-04 -1.536e-05  2.183e-04  1.750e-03 
#> 
#> Coefficients:
#>                             Estimate Std. Error   t value Pr(>|t|)    
#> (Intercept)                5.246e+00  3.027e-03  1732.815  < 2e-16 ***
#> distributed_lag_baro_0    -1.949e-01  1.203e-02   -16.201  < 2e-16 ***
#> distributed_lag_baro_1     2.365e-02  9.164e-03     2.581 0.009857 ** 
#> distributed_lag_baro_2     7.478e-03  5.490e-03     1.362 0.173165    
#> distributed_lag_baro_4     1.377e-02  2.583e-03     5.329 9.94e-08 ***
#> distributed_lag_baro_7     9.558e-03  1.279e-03     7.472 8.06e-14 ***
#> distributed_lag_baro_12    7.833e-03  5.735e-04    13.659  < 2e-16 ***
#> distributed_lag_baro_22    6.194e-03  2.442e-04    25.364  < 2e-16 ***
#> distributed_lag_baro_37    3.456e-03  1.009e-04    34.250  < 2e-16 ***
#> distributed_lag_baro_63    1.618e-03  4.268e-05    37.905  < 2e-16 ***
#> distributed_lag_baro_106   1.987e-04  1.882e-05    10.561  < 2e-16 ***
#> distributed_lag_baro_179   7.293e-05  8.501e-06     8.579  < 2e-16 ***
#> distributed_lag_baro_302  -7.000e-05  4.643e-06   -15.075  < 2e-16 ***
#> distributed_lag_baro_509  -7.968e-02  1.580e-03   -50.432  < 2e-16 ***
#> distributed_lag_baro_856   2.068e-01  4.099e-03    50.455  < 2e-16 ***
#> distributed_lag_baro_1440 -1.272e-01  2.519e-03   -50.499  < 2e-16 ***
#> earthtide_cos_1            5.839e-06  3.516e-06     1.661 0.096754 .  
#> earthtide_sin_1           -2.185e-06  3.551e-06    -0.615 0.538351    
#> earthtide_cos_2           -8.980e-05  3.624e-06   -24.784  < 2e-16 ***
#> earthtide_sin_2            3.769e-05  3.654e-06    10.315  < 2e-16 ***
#> earthtide_cos_3           -2.521e-05  3.069e-06    -8.215  < 2e-16 ***
#> earthtide_sin_3            3.247e-05  2.976e-06    10.911  < 2e-16 ***
#> earthtide_cos_4           -2.293e-04  8.237e-06   -27.833  < 2e-16 ***
#> earthtide_sin_4           -3.859e-05  7.148e-06    -5.399 6.74e-08 ***
#> earthtide_cos_5            3.306e-06  3.424e-06     0.966 0.334237    
#> earthtide_sin_5            1.281e-05  3.425e-06     3.740 0.000184 ***
#> earthtide_cos_6           -1.040e-05  3.866e-06    -2.690 0.007139 ** 
#> earthtide_sin_6            1.148e-05  3.881e-06     2.959 0.003083 ** 
#> earthtide_cos_7           -3.470e-06  2.159e-06    -1.607 0.107992    
#> earthtide_sin_7            1.585e-05  2.161e-06     7.336 2.24e-13 ***
#> earthtide_cos_8           -5.746e-05  2.712e-06   -21.185  < 2e-16 ***
#> earthtide_sin_8            4.885e-05  2.713e-06    18.006  < 2e-16 ***
#> earthtide_cos_9           -3.015e-04  2.809e-06  -107.347  < 2e-16 ***
#> earthtide_sin_9            2.906e-04  2.796e-06   103.947  < 2e-16 ***
#> earthtide_cos_10          -5.472e-06  3.208e-06    -1.706 0.088029 .  
#> earthtide_sin_10          -5.066e-05  3.191e-06   -15.878  < 2e-16 ***
#> earthtide_cos_11          -2.415e-04  6.753e-06   -35.766  < 2e-16 ***
#> earthtide_sin_11           8.373e-05  6.089e-06    13.749  < 2e-16 ***
#> earthtide_cos_12          -4.485e-07  2.640e-06    -0.170 0.865086    
#> earthtide_sin_12           3.140e-06  2.644e-06     1.187 0.235142    
#> datetime_ns_01            -2.450e-02  3.813e-05  -642.432  < 2e-16 ***
#> datetime_ns_02            -3.382e-02  5.298e-05  -638.423  < 2e-16 ***
#> datetime_ns_03            -4.747e-02  4.648e-05 -1021.392  < 2e-16 ***
#> datetime_ns_04            -6.418e-02  4.692e-05 -1367.631  < 2e-16 ***
#> datetime_ns_05            -7.759e-02  4.513e-05 -1719.314  < 2e-16 ***
#> datetime_ns_06            -8.960e-02  4.598e-05 -1948.683  < 2e-16 ***
#> datetime_ns_07            -1.035e-01  4.563e-05 -2268.022  < 2e-16 ***
#> datetime_ns_08            -1.185e-01  4.645e-05 -2550.455  < 2e-16 ***
#> datetime_ns_09            -1.299e-01  4.668e-05 -2781.825  < 2e-16 ***
#> datetime_ns_10            -1.392e-01  3.033e-05 -4588.477  < 2e-16 ***
#> datetime_ns_11            -1.621e-01  9.562e-05 -1695.121  < 2e-16 ***
#> datetime_ns_12            -1.513e-01  1.913e-05 -7911.386  < 2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error: 0.0003716 on 35229 degrees of freedom
#>   (1440 observations deleted due to missingness)
#> Multiple R-squared:  0.9999, Adjusted R-squared:  0.9999 
#> F-statistic: 9.426e+06 on 51 and 35229 DF,  p-value: < 2.2e-16
```
