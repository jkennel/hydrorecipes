#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Building Block for a Recipe --------------------------------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' R6 Class representing a step
#'
#' The `Step` class hold common info for each step.
#'
#' @inheritParams recipes::step_center
#' @inheritParams recipes::step_pca
#'
#' @param step_name the name of the step
#'
#'
#' @importFrom collapse fmean fsd fscale fsum fquantile fndistinct flag
#' @importFrom collapse missing_cases varying rowbind
#' @importFrom collapse qDF qM qF qTBL mctl
#' @importFrom earthtide calc_earthtide
#' @importFrom R6 R6Class
#' @importFrom Bessel BesselK BesselJ BesselI
#'
#' @export
Step <- R6Class(
  classname = "step",
  public = list(
    type = NULL, # check, add, remove, update/modify

    # base step
    terms = NULL,
    role = NULL,
    trained = FALSE,
    skip = FALSE,
    columns = NULL,
    step_name = NULL,
    keep_original_cols = TRUE,
    id = NULL,
    prefix = NULL,

    check = NULL,
    new_columns = c(),

    initialize = function(terms, ...) {
      if (!missing(terms)) {
        if (length(terms) == 1) {
          self$terms <- get_terms_and_symbols(c(terms))
        } else {
          self$terms <- get_terms_and_symbols(terms)
        }
      }

      dots <- c(...)
      self$role <- dots$role
      self$skip <- dots$skip
      self$keep_original_cols <- dots$keep_original_cols
      self$step_name <- dots$step_name
      self$type <- dots$type
      self$prefix <- dots$prefix


      # super specific values
      if (is.null(self$prefix)) {
        self$prefix <- gsub("step_", "", self$step_name)
      }
      self$id <- rand_id(self$prefix)

      invisible(self)
    },

    # these are the base methods - can be overwritten in individual steps
    prep = function(new_data, info) {
      nms <- names(new_data)
      self$columns <- get_terms_from_info(self$terms, nms, info)
      self$trained <- TRUE

      invisible(self)
    },
    bake = function() {
      invisible(self)
    },
    tidy = function(i) {

      # print(i)
      # print(self$columns)
      # print(self$new_columns)
      # print(self$role)

      if (is.null(self$new_columns)) {
        self$new_columns <- self$columns
      }

      data.frame(
        index       = i,
        variable    = self$columns,
        columns     = self$new_columns,
        role        = self$role,
        step_name   = self$step_name,
        id          = self$id,
        type        = self$type
      )
    },
    response = function(co) {
      n <- length(co)
      list(
        x = rep(NA_real_, n),
        variable = rep("coefficient", n),
        value = co,
        step_id = rep(self$id, n)
      )
    }

  )
)




# b <- 10
#
# a <- function(..., x = 'blah', y = 1, z = NULL) {
#   tmp <- c(`...` = enquos(...), environment())
#   print(get_env(environment()))
#   print(as.list(environment()))
# }
#
# l <- function(..., x = 'blah', y = 1, z = NULL) {
#   b
# }
#
# z <- a(b, x = 'ada')
#
# do.call(l, z)
