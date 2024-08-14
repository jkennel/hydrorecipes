
# Hsieh 1987

data('hsieh_1987_fig_2_3')

storativity <- 1e-07
transmissivity <- 1e-03
radius_well <- 0.05
frequency <- transmissivity / (hsieh_1987_fig_2_3[variable == 'gain' & S == storativity]$dimensionless_frequency * radius_well^2)
tau   <- 1 / frequency

hseih <- tidal_hsieh_1987(frequency, storativity, transmissivity, radius_well)
expect_equal(Mod(hseih$response),
             hsieh_1987_fig_2_3[variable == 'gain' & S == storativity]$response, tolerance = 0.004)


frequency <- transmissivity / (hsieh_1987_fig_2_3[variable=='phase' & S == storativity]$dimensionless_frequency * radius_well^2)
tau   <- 1 / frequency

hseih <- tidal_hsieh_1987(frequency, storativity, transmissivity, radius_well)
expect_equal(unwrap(Arg(hseih$response)) * 180 / pi,
             hsieh_1987_fig_2_3[variable == 'phase' & S == storativity]$response, tolerance = 0.01)


