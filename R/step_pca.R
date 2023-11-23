#' R6 Class
#'
#' `StepPca` adjust the central value to zero.
#' @inheritParams Step
#'
#' @export
StepPca <- R6Class(
  classname = 'step_pca',
  inherit = Step,

  public = list(
    pca_results = list(),
    n_comp = NA_integer_,
    na_rm = NA,
    center = NA,
    scale = NA,

    # step specific variables
    initialize = function(...,
                          role = "predictor",
                          skip = FALSE,
                          na_rm = TRUE,
                          n_comp = 3,
                          center = TRUE,
                          scale = TRUE,
                          keep_original_cols = FALSE) {

      # get function parameters to pass to parent
      step_name    <- "step_pca"
      type         <- 'modify'
      inputs <- c(
        as.list(rlang::quos(...)),
        rlang::env_get_list(env = environment(),
                            formalArgs(super$initialize)[-1])
      )
      do.call(super$initialize, inputs)

      self$na_rm <- na_rm
      self$n_comp <- n_comp
      self$center <- center
      self$scale <- scale

      invisible(self)
    },
    prep = function(new_data) {

      if(self$center) {
        self$center <- collapse::fmean(new_data,
                                       na.rm = self$na_rm)
      } else {
        self$center <- rep(0.0, length(new_data))
      }

      if(self$scale) {
        self$scale <- collapse::fsd(new_data,
                                    na.rm = self$na_rm)
      } else {
        self$scale <- rep(1.0, length(new_data))
      }

      self$pca_results <- pca_list_rotation_eigen(new_data,
                                                  center = self$center,
                                                  scale = self$scale,
                                                  n_comp = self$n_comp)
    },
    # subtract the central value from a column
    bake = function(new_data) {

      new_data <- collapse::qM(scale_list_param_eigen(new_data,
                                    center = self$center,
                                    scale = self$scale))
      new_data <- collapse::mctl(new_data %*% self$pca_results)

      width <- floor(log10(self$n_comp)) + 1

      names(new_data) <- paste0("PC", formatC(1:self$n_comp,
                                              width = width,
                                              format = "d",
                                              flag = "0"))
      new_data
    }
  )
)

