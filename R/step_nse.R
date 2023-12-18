#' R6 Class
#'
#' `StepNSE` adds variable vectors.
#' @param vars name of vars
#'
#' @inheritParams Step
#'
#' @export
StepNSE <- R6Class(
  classname = 'step_nse',
  inherit = Step,

  public = list(

    # step specific variables

    initialize = function(terms,
                          role = "predictor",
                          skip = FALSE,
                          keep_original_cols = FALSE,
                          ...) {

      # get function parameters to pass to parent
      step_name    <- "step_nse"
      type         <- 'modify'
      enq          <- NULL

      env_nms <- names(environment())
      super_nms <- formalArgs(super$initialize)
      env_nms_sub <- intersect(super_nms, env_nms)[-1L]
      print(env_nms_sub)
      # print(as.list(substitute(terms, environment())))

      # inputs <- c(
      #   as.list(rlang::quos(...)),
      #   rlang::env_get_list(env = environment(),
      #                       formalArgs(super$initialize)[-1L])
      # )
      # print(str(inputs))
      inputs2 <- c(
        as.list(substitute(terms, globalenv())),
        as.list(environment())[env_nms_sub]
      )
      print('here')
      do.call(super$initialize, inputs2)

      invisible(self)
    },
    bake = function(new_data) {
      return(list(x = rep(2.0, length(unlist(new_data)))))
    }

  )
)
