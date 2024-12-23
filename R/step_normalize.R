#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Adjust the dispersion and central value --------------------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
StepNormalize <- R6Class(
  classname = "step_normalize",
  inherit = Step,
  public = list(
    center = c(),
    scale = c(),
    na_rm = TRUE,
    # step specific variables
    initialize = function(terms,
                          role = "predictor",
                          na_rm = TRUE,
                          ...) {
      # get function parameters to pass to parent
      terms <- substitute(terms)
      env_list <- get_function_arguments()
      env_list$step_name <- "step_normalize"
      env_list$type <- "modify"
      super$initialize(
        terms = terms,
        env_list[names(env_list) != "terms"],
        ...
      )

      self$na_rm <- na_rm

      invisible(self)
    },
    prep = function(data) {


      # print("00000")
      # print(head(data))
      self$center <- collapse::fmean(data, na.rm = self$na_rm)

      self$scale <- collapse::fsd(data, na.rm = self$na_rm)
      self$scale <- 1.0 / self$scale

    },
    # subtract the central value from a column and divide by the standard deviation
    bake = function(s) {

      # print(head(s[["result"]][self$columns]))
      s[["result"]][self$columns] <- (s[["result"]][self$columns] %r-% self$center) %r*%
        (self$scale)

      return(NULL)
    }
  )
)
