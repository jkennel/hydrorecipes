eig  <- log_lags(15, 86401)
arma <- log_lags_arma(15, 86401)

expect_equivalent(as.numeric(eig),
                  as.numeric(arma),
                  info = "log_lags and log_lags_arma are equivalent")
