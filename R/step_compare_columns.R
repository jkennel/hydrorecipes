# This is experimental...The goal is to compare two columns and remove values
# that fall far from the expected value.
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

    prep = function(data) {

      self$column_values <- collapse::fsd(abs(diff(data[[2L]])),
                                          na.rm = self$na_rm)

    },

    bake = function(s) {


      self$result <- list(c(FALSE,
                            abs(diff(s[["result"]][[1L]])) > (self$column_values * self$n_sd)))

      self$result <- setNames(self$result,
                              paste0(self$id, "_",
                                     self$columns[1L], "_",
                                     self$columns[2]))

      return(NULL)

    }
  )
)
