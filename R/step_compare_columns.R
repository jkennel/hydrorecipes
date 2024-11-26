# This is experimental...The goal is to compare two columns and remove values
# that fall far from the estimated value.
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Filter based on column comparison --------------------------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
StepCompareColumns <- R6Class(
  classname = "step_compare_columns",
  inherit = Step,
  public = list(

    # step specific variables
    data = NULL,
    compare = NULL,
    n_sd = NULL,
    na_rm = NULL,
    column_values = NULL,

    initialize = function(data,
                          compare,
                          role = "predictor",
                          n_sd = 4,
                          na_rm = TRUE,
                          ...) {

      data <- deparse(substitute(data))
      compare <- deparse(substitute(compare))
      env_list <- get_function_arguments()
      env_list$step_name <- 'step_compare_columns'
      env_list$type <- 'add'
      super$initialize(terms = c(as.symbol(data),
                                 as.symbol(compare)),
                       env_list,
                       ...)

      self$n_sd <- n_sd
      self$data <- data
      self$na_rm <- na_rm
      self$compare <- compare
      self$columns <- c(data, compare)

      invisible(self)
    },

    prep = function(new_data, info) {

      super$prep(new_data, info)
      new_data <- unclass(new_data)[self$columns]
      self$column_values <- collapse::fsd(abs(diff(new_data[[2]])),
                                          na.rm = self$na_rm)

    },

    bake = function(new_data) {

      new_data <- unclass(new_data)[self$columns]
      ret <- list(c(FALSE, abs(diff(new_data[[1]])) > (self$column_values * self$n_sd)))

      self$result <- setNames(ret, paste0(self$id, "_", self$columns[1], "_", self$columns[2]))
      self$result

    }
  )
)
