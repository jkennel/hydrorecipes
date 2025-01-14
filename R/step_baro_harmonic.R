#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Calculate Barometric Efficiency using Harmonic tide methods ------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
StepBaroHarmonic <- R6Class(
  classname = "step_baro_rau",
  inherit = Step,
  public = list(

    water_level = NULL,
    barometric_pressure = NULL,
    earth_tide = NULL,
    time = NULL,

    frequency = NULL,
    cycle_size = NULL,
    start = NULL,
    inverse = NULL,

    barometric_efficiency = list(),

    # step specific variables
    initialize = function(time,
                          water_level,
                          barometric_pressure,
                          earth_tide,
                          frequency = c(1.9324, 2.0), # M2 and S2
                          cycle_size = 86400,
                          start = 0.0,
                          inverse = TRUE,
                          role = "augment",
                          ...) {
      # get function parameters to pass to parent
      time <- deparse(substitute(time))
      water_level <- deparse(substitute(water_level))
      barometric_pressure <- deparse(substitute(barometric_pressure))
      earth_tide <- deparse(substitute(earth_tide))
      env_list <- get_function_arguments()
      env_list$step_name <- "step_baro_harmonic"
      env_list$type <- "augment"
      super$initialize(
        terms = c(
          as.symbol(time),
          as.symbol(water_level),
          as.symbol(barometric_pressure),
          as.symbol(earth_tide)),
        env_list,
        ...
      )

      self$time <- time
      self$water_level <- water_level
      self$barometric_pressure <- barometric_pressure
      self$earth_tide <- earth_tide
      self$frequency <- frequency
      self$cycle_size <- cycle_size
      self$start <- start
      self$inverse <- inverse

      self$columns <- c(time, water_level, barometric_pressure, earth_tide)

      invisible(self)
    },
    bake = function(s) {

      # this is a hack to deal with NSE issues
      names(s[["result"]])[1L] <- "time_col"

      nms <- names(s[["result"]])

      # create regression formula
      formula_txt <- paste0(paste(nms[-1L], collapse = '+'), "~", nms[1])

      # include linear trend and intercept
      harmonics <- Recipe$new(formula = as.formula(formula_txt), s[["result"]])$
        add_step(StepIntercept$new())$
        add_step(StepHarmonic$new(time_col,
                                  frequency = self$frequency,
                                  cycle_size = self$cycle_size,
                                  starting_value = self$start))$
        plate("m")


      X <- harmonics[, -(2:4)]
      Y <- harmonics[, (2:4)]


      soln <- llt_solve(X, Y)

      co_names <- colnames(X)
      wh_s <- grep("sin", co_names)
      wh_c <- grep("cos", co_names)


      soln_cplx <- sin_cos_to_complex(c = soln[wh_c,], s = -soln[wh_s,])

      self$barometric_efficiency <- be_harmonic_cpp(soln_cplx, self$inverse)
      names(self$barometric_efficiency) <- c("ratio", "acworth", "rau")

      dt <- diff(as.numeric(new_data[[1L]])[1:2])
      cycle_size <- self$cycle_size / dt
      be_tf <- Mod(be_transfer(collapse::qM(new_data[2:4]),
                               5,
                               TRUE,
                               TRUE,
                               0.1,
                               2.0,
                               self$cycle_size / dt)[1L])

      names(be_tf) <- "tf"

      if (self$inverse) {
        be_tf[1L] <- 1.0 - be_tf[1L]
      }

      self$barometric_efficiency <- c(self$barometric_efficiency, be_tf)


      return(NULL)

    }
  )
)

sin_cos_to_complex <- function(c, s) {
  complex(real = t(c), imaginary = t(s))
}
