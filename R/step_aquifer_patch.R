#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Barker and Herbert Two Radial Patches Step -----------------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' R6 Class
#'
#' `StepAquiferPatch` barker_herbert 1982 solution for radial patches.
#'
#' @param time
#' @param flow_rate
#' @param thickness
#' @param radius
#' @param radius_patch
#' @param hydraulic_conductivity_inner
#' @param hydraulic_conductivity_outer
#' @param specific_storage_inner
#' @param specific_storage_outer
#'
#' @inheritParams Step
#'
#' @references
#' Barker, J.A., and R. Herbert, 1982: Pumping tests in patchy
#'  aquifers, Ground Water, vol. 20, No. 2, pp. 150-155.
#'
#' Butler, J.J., 1988: Pumping tests in nonuniform aquifers – The radially
#'  symmetric case, Journal of Hydrology, Vol. 101, pp. 15-30.
#'
#' @family aquifer
#'
#' @export
StepAquiferPatch <- R6Class(
  classname = "step_aquifer_patch",
  inherit = Step,
  public = list(

    # step specific variables
    time = NULL,
    flow_rate = NULL,
    thickness = NULL,
    radius = NULL,
    radius_patch = NULL,
    specific_storage_inner = NULL,
    specific_storage_outer = NULL,
    hydraulic_conductivity_inner = NULL,
    hydraulic_conductivity_outer = NULL,

    # calculated
    storativity_inner = NULL,
    storativity_outer = NULL,
    transmissivity_inner = NULL,
    transmissivity_outer = NULL,
    n_terms = NULL,
    initialize = function(time,
                          flow_rate = 0.01,
                          thickness = 1.0,
                          radius = 200.0,
                          radius_patch = 100.0,
                          specific_storage_inner = 1.0e-6,
                          specific_storage_outer = 1.0e-5,
                          hydraulic_conductivity_inner = 1.0e-4,
                          hydraulic_conductivity_outer = 1.0e-6,
                          n_terms = 12L,
                          role = "predictor",
                          ...) {
      # get function parameters to pass to parent
      time <- deparse(substitute(time))
      env_list <- get_function_arguments()
      env_list$step_name <- "step_aquifer_patch"
      env_list$type <- "add"
      super$initialize(
        terms = c(as.symbol(time), as.symbol(flow_rate)),
        env_list
      )

      # step specific values
      self$time <- time
      self$columns <- self$time
      self$flow_rate <- flow_rate

      self$thickness <- thickness
      self$radius <- radius
      self$radius_patch <- radius_patch

      # K & Ss
      self$specific_storage_inner <- specific_storage_inner
      self$specific_storage_outer <- specific_storage_outer
      self$hydraulic_conductivity_inner <- hydraulic_conductivity_inner
      self$hydraulic_conductivity_outer <- hydraulic_conductivity_outer

      # T & S
      self$storativity_inner <- specific_storage_inner * thickness
      self$storativity_outer <- specific_storage_outer * thickness
      self$transmissivity_inner <- hydraulic_conductivity_inner * thickness
      self$transmissivity_outer <- hydraulic_conductivity_outer * thickness

      self$n_terms <- n_terms


      invisible(self)
    },
    # (frecipes:::stehfest_barker_herbert(as.numeric(1:n),1.0, 100.0, 200.0, 1e-3, 1e-3, 1e-5, 1e-5, 12L)[[1]])

    bake = function(new_data) {


      # Eigen::VectorXd time,
      # double radius,
      # double radius_patch,
      # double t_1,
      # double t_2,
      # double s_1,
      # double s_2,
      # double Q,
      # double prec,
      # int n_terms

      bh <- setNames(list(barker_herbert(
        time = new_data[[1]],
        radius = self$radius,
        radius_patch = self$radius_patch,
        t_1 = self$transmissivity_inner,
        t_2 = self$transmissivity_outer,
        s_1 = self$storativity_inner,
        s_2 = self$storativity_outer,
        Q = self$flow_rate,
        1e-8,
        n_terms = self$n_terms
        )), self$id)

      return(bh)


      # setNames(
      #   gwr_barker_herbert(
      #     new_data[[1]],
      #     self$flow_rate,
      #     self$radius,
      #     self$radius_patch,
      #     self$transmissivity_inner,
      #     self$transmissivity_outer,
      #     self$storativity_inner,
      #     self$storativity_outer,
      #     self$n_terms
      #   ),
      #   self$id
      # )
    }
  )
)
