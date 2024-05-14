#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Divide a Term into Intervals and do Dummy Encoding ---------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
StepFindInterval <- R6Class(
  classname = "step_find_interval",
  inherit = Step,
  public = list(

    # step specific variables
    vec = NULL,
    n_vec = NULL,

    initialize = function(terms,
                          vec,
                          role = "predictor",
                          ...) {
      # get function parameters to pass to parent
      terms <- substitute(terms)
      env_list <- get_function_arguments()
      env_list$step_name <- "step_find_interval"
      env_list$type <- "add"
      super$initialize(
        terms = terms,
        env_list[names(env_list) != "terms"],
        ...
      )



      # step specific values
      self$vec <- sort(vec)
      self$n_vec <- length(vec)

      invisible(self)
    },
    bake = function(new_data) {
      column_name <- self$columns

      self$new_columns <- c()

      dum <- list()
      for (i in seq_along(column_name)) {
        dum[[i]] <- to_dummy_list(unclass(new_data)[[i]], self$vec)

        nn <- name_columns(
          self$prefix,
          column_name[i],
          length(dum[[i]])
        )

        names(dum[[i]]) <- nn
        self$new_columns <- c(self$new_columns, nn)
      }
      unlist(dum, recursive = FALSE)
    }
  )
)
