#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Convert a Column to Dummy Encoding Step --------------------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
StepDummy <- R6Class(
  classname = "step_dummy",
  inherit = Step,
  public = list(

    # step specific variables
    levels = NULL,
    one_hot = NULL,

    initialize = function(terms,
                          one_hot = FALSE,
                          role = "predictor",
                          ...) {

      # get function parameters to pass to parent
      terms <- substitute(terms)
      env_list <- get_function_arguments()
      env_list$step_name <- "step_dummy"
      env_list$type <- "add"
      super$initialize(
        terms = terms,
        env_list[names(env_list) != "terms"],
        ...
      )

      invisible(self)

      self$one_hot <- one_hot

    },
    prep = function(new_data, info) {

      super$prep(new_data, info)
      self$levels <- lapply(unclass(new_data)[self$columns], levels)
      invisible(self)

    },
    bake = function(new_data) {
      column_name <- self$columns

      dum <- list()
      for (i in seq_along(column_name)) {
        # check for new level(s)
        if (sum(levels(unclass(new_data)[[i]]) %!in% self$levels[[i]]) > 0L) {
          warning(paste0("New levels found during bake step. (", self$id, ")"))
        }

        dum[[i]] <- to_dummy(unclass(new_data)[[i]], self$one_hot)

        names(dum[[i]]) <- name_columns(
          self$id,
          column_name[i],
          n = length(dum[[i]])
        )
      }

      self$result <- unlist(dum, recursive = FALSE)
      self$result

    }
  )
)
