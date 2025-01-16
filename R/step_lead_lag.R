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
      self$lag <- as.integer(lag[order(lag)])
      self$n_shift <- as.integer(n_shift)
      self$n_subset <- as.integer(n_subset)

      invisible(self)
    },
    bake = function(s) {

      ll <- list()

      self$new_columns <- c()
      for (i in seq_along(self$columns)) {
        column_name <- self$columns[i]
        # if (self$n_subset == 1) {
        #   ll[[i]] <- collapse::flag(s[["result"]][column_name], self$lag)
        # } else {
          ll[[i]] <- lag_list(s[["result"]][[column_name]],
                              self$lag,
                              n_subset = self$n_subset,
                              n_shift = self$n_shift
          )
        # }

        # print(str(ll))
        nn <- name_columns(
          self$prefix,
          column_name,
          n = length(self$lag)
        )

        names(ll[[i]]) <- nn
        self$new_columns <- c(self$new_columns, nn)
      }

      self$result <- unlist(ll, recursive = FALSE)
      return(NULL)

    },
    response = function(co) {


      nr <- nrow(co)
      nc <- ncol(co)


      list(x = rep.int(self$lag, 2L * nc),
           variable = rep(c("coefficient", "cumulative"), each = nr * nc),
           value = c(co, collapse::fcumsum(co)),
           step_id = rep.int(self$id, 2L * nr * nc),
           outcome = rep(rep(colnames(co), each = nr), 2L),
           term = rep.int(self$prefix, 2L * nr * nc),
           step_columns = rep.int(paste(self$columns, collapse = "_"), 2L * nr * nc))




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
