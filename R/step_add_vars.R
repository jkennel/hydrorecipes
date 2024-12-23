#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Add variables to Recipe Step -------------------------------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
StepAddVars <- R6Class(

  classname = "step_add_vars",
  inherit = Step,

  public = list(

    # step specific variables
    initialize = function(terms,
                          role = "predictor",
                          ...) {

      # get function parameters to pass to parent - strings can be passed
      is_string <- tryCatch(as.character(terms),
                            error = function(terms) (terms))


      if ("error" %in% class(is_string)) {
        terms <- (substitute(terms))
      } else if (!class(terms) == "character") {
        terms <- substitute(terms)
      }

      env_list <- get_function_arguments()
      env_list$step_name <- "step_add_vars"
      env_list$type <- "add_from_template"

      super$initialize(
        terms = terms,
        env_list[names(env_list) != "terms"],
        ...
      )

      self$new_columns <- as.character(terms)


      invisible(self)
    }

  )
)

