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
    # aquifer thickness
    thickness = NULL,
    # leakage hantush leakage
    leakage = NULL,
    # radius distance to monitoring interval
    radius = NULL,
    # specific_storage aquifer specific_storage
    specific_storage = NULL,
    # hydraulic_conductivity aquifer hydraulic_conductivity
    hydraulic_conductivity = NULL,
    # max_terms number of terms to use in Hantush solution.  More is more
    #   precise but slower.
    precision = NULL,
    initialize = function(time,
                          flow_rate,
                          thickness = 1.0,
                          leakage = 100.0,
                          radius = 100.0,
                          specific_storage = 1e-6,
                          hydraulic_conductivity = 1e-4,
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
      self$thickness <- thickness
      self$specific_storage <- specific_storage
      self$hydraulic_conductivity <- hydraulic_conductivity
      self$precision <- precision

      self$columns <- c(time, flow_rate)

      invisible(self)
    },
    bake = function(s) {

      self$new_columns <- self$prefix
      # self$columns <- paste(self$columns, collapse = ",")

      hj <- hantush_jacob(
        s[["result"]][[self$columns[[1L]]]],
        s[["result"]][[self$columns[[2L]]]],
        self$radius,
        self$specific_storage * self$thickness,
        self$hydraulic_conductivity * self$thickness,
        self$leakage,
        self$precision
      )


      self$result <- setNames(hj, self$new_columns)
      return(NULL)
    }
  )
)
