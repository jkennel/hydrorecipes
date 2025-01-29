#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Cooper-Bredehoeft-Papadopulos Slug Test --------------------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
StepSlugCbp <- R6Class(
  classname = "step_aquifer_constant_drawdown",
  inherit = Step,
  public = list(

    # step specific variables
    time = NULL,
    thickness = NULL,
    radius = NULL,
    radius_casing = NULL,
    radius_well = NULL,
    specific_storage = NULL,
    hydraulic_conductivity = NULL,
    head_0 = NULL,
    n_terms = NULL,

    initialize = function(time,
                          radius = 0.1,
                          radius_casing = 0.1,
                          radius_well = 0.1,
                          specific_storage = 1.0e-6,
                          hydraulic_conductivity = 1.0e-4,
                          head_0 = 1.0,
                          thickness = 1.0,
                          n_terms = 16,
                          role = "predictor",
                          ...) {

      # get function parameters to pass to parent
      time <- deparse(substitute(time))
      env_list <- get_function_arguments()
      env_list$step_name <- "step_slug_cbp"
      env_list$type <- "add"
      super$initialize(
        terms = as.symbol(time),
        env_list,
        ...
      )


      # step specific values
      self$radius <- radius
      self$radius_casing <- radius_casing
      self$radius_well <- radius_well
      self$specific_storage <- specific_storage
      self$hydraulic_conductivity <- hydraulic_conductivity
      self$head_0 <- head_0
      self$thickness <- thickness

      self$n_terms <- n_terms
      self$columns <- time

      invisible(self)
    },
    bake = function(s) {

      Tr  <- self$hydraulic_conductivity * self$thickness
      S  <- self$specific_storage * self$thickness

      # print(self$radius)
      # print(self$radius_casing)
      # print(self$radius_well)
      # print(self$hydraulic_conductivity)
      # print(self$specific_storage)
      # print(self$thickness)
      # print(self$head_0)
      # print(self$n_terms)
      cbp <- list(cooper_bredehoeft_papadopulos_laplace(
        time = s[["result"]][[self$columns]],
        r = self$radius,
        r_c = self$radius_casing,
        r_w = self$radius_well,
        Tr = Tr,
        S = S,
        h_0 = self$head_0,
        # prec = 1e-8,
        n_terms = self$n_terms))

      self$new_columns <- self$prefix

      self$result <- setNames(cbp, self$new_columns)

      return(NULL)
    }
  )
)
