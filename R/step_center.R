#' R6 Class
#'
#' `StepCenter` adjust the central value to zero.
#' @inheritParams Step
#'
#' @export
StepCenter <- R6Class(
  classname = 'step_center',
  inherit = Step,


  public = list(
    means = c(),
    na_rm = NA,
    center_fun = NULL,

    # step specific variables
    initialize = function(...,
                          role = "modify",
                          skip = FALSE,
                          na_rm = TRUE,
                          center_fun = collapse::fmean,
                          keep_original_cols = FALSE) {

      # get function parameters to pass to parent
      step_name    <- "step_center"
      type         <- 'modify'
      inputs <- c(
        as.list(rlang::quos(...)),
        rlang::env_get_list(env = environment(),
                            formalArgs(super$initialize)[-1])
      )
      do.call(super$initialize, inputs)

      self$na_rm <- na_rm
      self$center_fun <- center_fun

      invisible(self)
    },
    prep = function(new_data) {
      self$means <- vapply(new_data[self$columns],
                           FUN = self$center_fun,
                           FUN.VALUE = numeric(1),
                           na.rm = self$na_rm)
    },
    bake = function(new_data) {

      for(i in seq_along(self$columns)) {
        new_data[self$columns[i]] <- NULL
      }

      new_data
    }
  )
)

