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
#' @importFrom rlang quos enquos env_get_list
#' @importFrom collapse fmean fsd fscale fsum fquantile fndistinct flag
#' @importFrom collapse missing_cases
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

      # super specific values
      self$prefix <- gsub("step_", "", self$step_name)
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
    tidy = function() {
      if (self$type == "add") {
        data.frame(
          step_name     = self$step_name,
          id            = self$id,
          columns       = rep(self$columns, each = sapply(self$result, ncol)),
          columns_added = self$new_columns,
          type          = self$type,
          role          = self$role
        )
      } else {

      }
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
