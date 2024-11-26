#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Create lagged or leaded terms ------------------------------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
StepLeadLag <- R6Class(
  classname = "step_lead_lag",
  inherit = Step,
  public = list(

    # step specific variables
    lag = NULL,
    n_shift = NULL,
    n_subset = NULL,
    initialize = function(terms,
                          lag,
                          n_shift = 0L,
                          n_subset = 1L,
                          role = "predictor",
                          ...) {
      # get function parameters to pass to parent
      terms <- substitute(terms)
      env_list <- get_function_arguments()
      env_list$step_name <- "step_lead_lag"
      env_list$type <- "add"
      super$initialize(
        terms = terms,
        env_list[names(env_list) != "terms"],
        ...
      )


      # step specific values
      self$lag <- as.integer(sort(lag))
      self$n_shift <- as.integer(n_shift)
      self$n_subset <- as.integer(n_subset)

      invisible(self)
    },
    bake = function(new_data) {
      column_name <- self$columns

      self$new_columns <- c()

      ll <- list()
      for (i in seq_along(column_name)) {
        if (self$n_subset == 1) {
          ll[[i]] <- collapse::flag(new_data[i], self$lag)
        } else {
          ll[[i]] <- lag_list(unclass(new_data)[[i]],
                              self$lag,
                              n_subset = self$n_subset,
                              n_shift = self$n_shift
          )
        }

        nn <- name_columns(
          self$prefix,
          column_name[i],
          n = length(self$lag)
        )

        names(ll[[i]]) <- nn
        self$new_columns <- c(self$new_columns, nn)
      }

      self$result <- unlist(ll, recursive = FALSE)
      self$result

    },
    response = function(co) {


      nr <- nrow(co)
      nc <- ncol(co)


      list(x = rep(self$lag, 2L * nc),
           variable = rep(c("coefficient", "cumulative"), each = nr * nc),
           value = c(co, collapse::fcumsum(co)),
           step_id = rep(self$id, 2L * nr * nc),
           outcome = rep(rep(colnames(co), each = nr), 2L),
           term = rep("lead_lag", 2L * nr * nc))




      # n <- length(co)
      #
      # variable <- c(
      #   rep("coefficient", n),
      #   rep("cumulative", n)
      # )
      #
      # value <- c(
      #   co,
      #   collapse::fcumsum(co)
      # )
      #
      # list(x = rep(self$lag, 2L),
      #      variable = variable,
      #      value = value,
      #      step_id = rep(self$id, 2L * n),
      #      term = "lead_lag")
    }
  )
)
