
<!-- README.md is generated from README.Rmd. Please edit that file -->

# hydrorecipes

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/hydrorecipes)](https://CRAN.R-project.org/package=hydrorecipes)
[![R-CMD-check](https://github.com/jkennel/hydrorecipes/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jkennel/hydrorecipes/actions/workflows/R-CMD-check.yaml)
[![License: GPL
v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

The goal of **hydrorecipes** is to supplement the [recipes
package](https://recipes.tidymodels.org) with a few steps that can help
deal with moderately sized water level datasets. These were developed
primarily with regression convolution in mind to be used with *lm* or
*glmnet*. The following steps are currently available:

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

You can install the development version of **hydrorecipes** from
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
a distributed lag model with 15 regressor terms. The *knots* are
logarithmically separated over two days to accurately capture early and
late time responses which can be caused by different physical
mechanisms.

``` r
library(hydrorecipes)
library(earthtide)
library(tidyr)
library(ggplot2)

data(transducer)

# convert to numeric because step_ns doesn't handle POSIXct
transducer$datetime <- as.numeric(transducer$datetime)

unique(diff(transducer$datetime)) # times are regularly spaced
#> [1] 120

# Earth tide inputs
wave_groups <- earthtide::eterna_wavegroups
wave_groups <- na.omit(wave_groups[wave_groups$time == '1 month', ])
wave_groups <- wave_groups[wave_groups$start > 0.5,]
latitude  <- 34.0
longitude <- -118.5
elevation <- 500

# create recipe 
rec <- recipe(wl~baro+datetime, transducer) |>
  step_distributed_lag(baro, knots = log_lags(15, 86400 * 2 / 120)) |>
  step_earthtide(datetime,
                 latitude = latitude,
                 longitude = longitude,
                 elevation = elevation,
                 astro_update = 1,
                 wave_groups = wave_groups) |>
  step_ns(datetime, deg_free = 15) |>
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
#> -1.704e-03 -2.225e-04 -1.595e-05  2.001e-04  1.434e-03 
#> 
#> Coefficients:
#>                             Estimate Std. Error   t value Pr(>|t|)    
#> (Intercept)                5.193e+00  2.901e-03  1790.116  < 2e-16 ***
#> distributed_lag_baro_0    -1.947e-01  1.065e-02   -18.276  < 2e-16 ***
#> distributed_lag_baro_1     2.404e-02  8.116e-03     2.962 0.003061 ** 
#> distributed_lag_baro_2     7.380e-03  4.861e-03     1.518 0.129010    
#> distributed_lag_baro_4     1.391e-02  2.288e-03     6.078 1.23e-09 ***
#> distributed_lag_baro_7     9.552e-03  1.133e-03     8.433  < 2e-16 ***
#> distributed_lag_baro_12    7.935e-03  5.078e-04    15.625  < 2e-16 ***
#> distributed_lag_baro_22    6.132e-03  2.163e-04    28.354  < 2e-16 ***
#> distributed_lag_baro_37    3.510e-03  8.935e-05    39.278  < 2e-16 ***
#> distributed_lag_baro_63    1.592e-03  3.779e-05    42.143  < 2e-16 ***
#> distributed_lag_baro_106   2.312e-04  1.665e-05    13.886  < 2e-16 ***
#> distributed_lag_baro_179   5.677e-05  7.512e-06     7.558 4.20e-14 ***
#> distributed_lag_baro_302  -4.676e-05  4.107e-06   -11.385  < 2e-16 ***
#> distributed_lag_baro_509  -7.841e-02  1.399e-03   -56.033  < 2e-16 ***
#> distributed_lag_baro_856   2.035e-01  3.630e-03    56.046  < 2e-16 ***
#> distributed_lag_baro_1440 -1.251e-01  2.231e-03   -56.078  < 2e-16 ***
#> earthtide_cos_1           -1.672e-05  3.104e-06    -5.386 7.24e-08 ***
#> earthtide_sin_1           -2.621e-06  3.135e-06    -0.836 0.403216    
#> earthtide_cos_2            8.041e-05  3.211e-06    25.041  < 2e-16 ***
#> earthtide_sin_2           -4.756e-05  3.239e-06   -14.684  < 2e-16 ***
#> earthtide_cos_3            1.838e-05  2.747e-06     6.692 2.23e-11 ***
#> earthtide_sin_3           -2.210e-05  2.683e-06    -8.239  < 2e-16 ***
#> earthtide_cos_4            2.321e-04  7.297e-06    31.805  < 2e-16 ***
#> earthtide_sin_4            1.630e-05  6.347e-06     2.568 0.010226 *  
#> earthtide_cos_5           -6.898e-06  3.003e-06    -2.297 0.021617 *  
#> earthtide_sin_5           -1.605e-05  3.004e-06    -5.341 9.32e-08 ***
#> earthtide_cos_6            1.615e-05  3.432e-06     4.705 2.55e-06 ***
#> earthtide_sin_6           -1.245e-05  3.445e-06    -3.615 0.000301 ***
#> earthtide_cos_7            3.086e-06  1.907e-06     1.618 0.105654    
#> earthtide_sin_7           -1.594e-05  1.909e-06    -8.348  < 2e-16 ***
#> earthtide_cos_8            5.805e-05  2.366e-06    24.535  < 2e-16 ***
#> earthtide_sin_8           -4.849e-05  2.367e-06   -20.488  < 2e-16 ***
#> earthtide_cos_9            3.043e-04  2.490e-06   122.224  < 2e-16 ***
#> earthtide_sin_9           -2.929e-04  2.478e-06  -118.172  < 2e-16 ***
#> earthtide_cos_10          -8.925e-06  3.062e-06    -2.915 0.003565 ** 
#> earthtide_sin_10           5.301e-05  3.054e-06    17.357  < 2e-16 ***
#> earthtide_cos_11           2.210e-04  5.982e-06    36.938  < 2e-16 ***
#> earthtide_sin_11          -7.191e-05  5.399e-06   -13.319  < 2e-16 ***
#> earthtide_cos_12           1.078e-06  2.338e-06     0.461 0.644709    
#> earthtide_sin_12          -1.086e-06  2.342e-06    -0.464 0.642922    
#> datetime_ns_01            -2.149e-02  5.417e-05  -396.716  < 2e-16 ***
#> datetime_ns_02            -3.292e-02  6.895e-05  -477.473  < 2e-16 ***
#> datetime_ns_03            -3.858e-02  6.018e-05  -640.965  < 2e-16 ***
#> datetime_ns_04            -5.205e-02  7.135e-05  -729.510  < 2e-16 ***
#> datetime_ns_05            -6.417e-02  6.298e-05 -1018.779  < 2e-16 ***
#> datetime_ns_06            -7.597e-02  6.574e-05 -1155.742  < 2e-16 ***
#> datetime_ns_07            -8.600e-02  6.383e-05 -1347.281  < 2e-16 ***
#> datetime_ns_08            -9.567e-02  6.486e-05 -1475.032  < 2e-16 ***
#> datetime_ns_09            -1.074e-01  6.477e-05 -1658.190  < 2e-16 ***
#> datetime_ns_10            -1.188e-01  6.598e-05 -1800.144  < 2e-16 ***
#> datetime_ns_11            -1.294e-01  6.432e-05 -2011.838  < 2e-16 ***
#> datetime_ns_12            -1.376e-01  6.325e-05 -2175.965  < 2e-16 ***
#> datetime_ns_13            -1.453e-01  4.279e-05 -3395.285  < 2e-16 ***
#> datetime_ns_14            -1.692e-01  1.362e-04 -1242.976  < 2e-16 ***
#> datetime_ns_15            -1.532e-01  2.308e-05 -6640.503  < 2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error: 0.0003291 on 35226 degrees of freedom
#>   (1440 observations deleted due to missingness)
#> Multiple R-squared:  0.9999, Adjusted R-squared:  0.9999 
#> F-statistic: 1.135e+07 on 54 and 35226 DF,  p-value: < 2.2e-16
```

## Decomposition

``` r
pred <- predict_terms(fit = fit, 
                      rec = rec,
                      data = input)
pred <- bind_cols(transducer[, c('datetime', 'wl')], pred)
pred_long <- pivot_longer(pred, cols = !datetime)
ggplot(pred_long, aes(x = datetime, y = value)) +
  geom_line() + 
  facet_grid(name~., scales = 'free_y') + 
  theme_bw()
#> Warning: Removed 1440 row(s) containing missing values (geom_path).
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="90%" />

## Response

``` r
resp <- response(fit, rec)
resp_ba <- resp[resp$name == 'cumulative',]
resp_ba <- resp_ba[resp_ba$term == 'baro',]

ggplot(resp_ba, aes(x = x * 120, y = value)) +
  ggtitle('Barometric Loading response') + 
  xlab('lag (seconds)') +
  ylab('Cumulative response') +
  scale_y_continuous(limits = c(0,1)) +
  geom_line() + 
  theme_bw()
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="90%" />

``` r

resp_et <- resp[resp$name %in% c('amplitude', 'phase'),]
ggplot(resp_et, aes(x = x, xend = x, y = 0, yend = value)) +
  geom_segment() + 
  ggtitle('Earthtide Response') +
  xlab('Frequency (cycles per day)') +
  facet_grid(name~., scales = 'free_y') + 
  theme_bw()
```

<img src="man/figures/README-unnamed-chunk-5-2.png" width="90%" />
