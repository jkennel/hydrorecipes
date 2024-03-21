kern_conf <- 0.38
max_syn_lag <- 86400

lags <- list(time = as.numeric(0:max_syn_lag))


# vadose kernels ----------------------------------------------------------
kern_vad_a <- recipe(time~., data = unclass(lags)) |>
  step_vadose_weeks(time = time,
                    air_diffusivity = 0.1,
                    thickness = 40,
                    precision = 1e-16,
                    inverse = TRUE) |>
  plate()
