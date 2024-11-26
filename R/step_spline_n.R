#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# N-Spline Step ----------------------------------------------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
StepSplineN <- R6Class(
  classname = "step_spline_n",
  inherit = Step,
  public = list(

    # step specific variables
    df = 0L,
    internal_knots = NULL,
    degree = 3L,
    intercept = FALSE,
    boundary_knots = NULL,
    periodic = FALSE,
    initialize = function(terms,
                          df = 0L,
                          internal_knots = NULL,
                          boundary_knots = NULL,
                          intercept = FALSE,
                          periodic = FALSE,
                          degree = 3L,
                          role = "predictor",
                          ...) {
      # get function parameters to pass to parent
      terms <- substitute(terms)
      env_list <- get_function_arguments()
      env_list$step_name <- "step_spline_n"
      env_list$type <- "add"
      super$initialize(
        terms = terms,
        env_list[names(env_list) != "terms"]
      )

      # step specific values
      if (df != 0L) {
        self$df <- df + (1L - intercept) + 2L
      } else {
        self$df <- 0L
      }

      self$internal_knots <- as.numeric(sort(internal_knots))
      self$degree <- degree
      self$intercept <- intercept
      self$boundary_knots <- as.numeric(boundary_knots)
      self$periodic <- periodic

      invisible(self)
    },
    prep = function(new_data, info) {
      super$prep(new_data, info)

      if (self$df != 0L) {
        ik <- collapse::fquantile(unclass(new_data)[[self$columns]],
                                  probs = seq(0.0, 1.0, length.out = self$df - 2),
                                  na.rm = TRUE
        )
        self$boundary_knots <- ik[c(1L, length(ik))]
        self$internal_knots <- ik[-c(1L, length(ik))]
      }
    },
    bake = function(new_data) {
      column_name <- self$columns

      self$new_columns <- c()

      basis <- list()
      for (i in seq_along(column_name)) {
        basis[[i]] <- n_spline_list(
          x = unclass(new_data)[[i]],
          df = self$df,
          degree = 3L,
          internal_knots = self$internal_knots,
          boundary_knots = self$boundary_knots,
          complete_basis = self$intercept
        )

        nn <- name_columns(
          self$prefix,
          column_name,
          n = length(basis[[i]])
        )
        names(basis[[i]]) <- nn
        self$new_columns <- c(self$new_columns, nn)
      }

      self$result <- unlist(basis, recursive = FALSE)
      self$result

    }
  )
)
