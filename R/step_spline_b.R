#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# B-Spline Step ----------------------------------------------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
StepSplineB <- R6Class(
  classname = "step_spline_b",
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
      env_list$step_name <- "step_spline_b"
      env_list$type <- "add"
      super$initialize(
        terms = terms,
        env_list[names(env_list) != "terms"]
      )

      # step specific values
      if (df != 0L) {
        self$df <- df + (1L - intercept)
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
    prep = function(data) {

      if (self$df != 0L) {
        ik <- collapse::fquantile(data[[1L]],
                                  probs = seq(0.0, 1.0, length.out = self$df - 2L),
                                  na.rm = TRUE
        )
        self$boundary_knots <- ik[c(1L, length(ik))]
        self$internal_knots <- ik[-c(1L, length(ik))]
      }
    },
    bake = function(s) {

      self$new_columns <- c()

      basis <- list()
      for (i in seq_along(self$columns)) {
        column_name <- self$columns[i]

        basis[[i]] <- b_spline_list(
          x = s[["result"]][[column_name]],
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

      return(NULL)
    }
  )
)
