#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# FFT Output divided by kernel -------------------------------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
StepKernelDivideNaive <- R6Class(
  classname = "step_kernel_divide",
  inherit = Step,
  public = list(

    # step specific variables
    kernel = NULL,
    initialize = function(terms,
                          kernel,
                          role = "predictor",
                          ...) {

      # get function parameters to pass to parent
      terms <- substitute(terms)
      env_list <- get_function_arguments()
      env_list$step_name <- "step_kernel_filter_naive"
      env_list$type <- "add"
      super$initialize(
        terms = terms,
        env_list[names(env_list) != "terms"],
        ...
      )


      # step specific values
      self$kernel <- if (!inherits(kernel, "list")) list(kernel) else kernel
      n_kernel <- length(kernel)

      invisible(self)
    },
    bake = function(s) {

      self$new_columns <- c()
      filt <- list()

      for (i in seq_along(self$columns)) {

        column_name <- self$columns[i]

        filt[[i]] <- convolve_divide_naive_list(
          s[["result"]][[column_name]],
          self$kernel
        )

        nn <-  name_columns(
          self$prefix,
          column_name,
          n = length(self$kernel)
        )

        names(filt[[i]]) <- nn
        self$new_columns <- c(self$new_columns, nn)

      }

      self$result <- unlist(filt, recursive = FALSE)

      return(NULL)
      # self$result
    }
  )
)
