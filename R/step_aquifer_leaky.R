#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Hantush and Jacob Leaky Aquifer Step -----------------------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
StepAquiferLeaky <- R6Class(
  classname = "step_aquifer_leaky",
  inherit = Step,
  public = list(

    # step specific variables
    # time prediction times
    time = NULL,
    # flow_rate well flow rates
    flow_rate = NULL,
    # leakage hantush leakage
    leakage = NULL,
    # radius distance to monitoring interval
    radius = NULL,
    # storativity aquifer storativity
    storativity = NULL,
    # transmissivity aquifer transmissivity
    transmissivity = NULL,
    # max_terms number of terms to use in Hantush solution.  More is more
    #   precise but slower.
    precision = NULL,
    initialize = function(time,
                          flow_rate,
                          leakage = 100.0,
                          radius = 100.0,
                          storativity = 1e-6,
                          transmissivity = 1e-4,
                          precision = 1e-10,
                          role = "predictor",
                          ...) {
      # get function parameters to pass to parent
      time <- deparse(substitute(time))
      flow_rate <- deparse(substitute(flow_rate))
      env_list <- get_function_arguments()
      env_list$step_name <- "step_aquifer_leaky"
      env_list$type <- "add"
      super$initialize(
        terms = c(as.symbol(time), as.symbol(flow_rate)),
        env_list,
        ...
      )


      # step specific values
      self$time <- time
      self$flow_rate <- flow_rate

      self$leakage <- leakage
      self$radius <- radius

      # K & Ss
      self$storativity <- storativity
      self$transmissivity <- transmissivity
      self$precision <- precision

      self$columns <- c(time, flow_rate)

      invisible(self)
    },
    bake = function(s) {

      self$new_columns <- self$prefix
      self$columns <- paste(self$columns, collapse = ",")

      hj <- hantush_jacob(
        s[["result"]][[self$columns[[1L]]]],
        s[["result"]][[self$columns[[2L]]]],
        self$radius,
        self$storativity,
        self$transmissivity,
        self$leakage,
        self$precision
      )


      self$result <- setNames(hj, self$new_columns)
      return(NULL)
    }
  )
)
