#' R6 Class representing a step
#'
#' The `Step` class hold common info for each step.
#'
#' @param type
#' @param role
#' @param trained
#' @param skip
#' @param columns
#' @param names_vector
#' @param new_columns
#' @param step_name
#' @param keep_original_cols
#' @param id
#' @param prefix
#' @param result
#'
#' @export
Step <- R6Class(

  classname = 'step',

  public = list(
    type = NULL,  # check, add, remove, update/modify

    # base step
    terms = NULL,
    role = NULL,
    trained = FALSE,
    skip = FALSE,
    columns = NULL,
    names_vector = NULL,
    new_columns = NULL,
    step_name = NULL,
    keep_original_cols = TRUE,
    id = NULL,
    prefix = NULL,

    # result should be in recipe
    result = NULL,

    initialize = function(..., role, skip, type,
                          keep_original_cols, step_name) {
      # super specific values
      self$terms   <- enquos(...)
      self$columns <- get_terms(self$terms)
      self$role    <- role
      self$skip    <- skip
      self$prefix  <- gsub("step_", "", step_name)
      self$id      <- rand_id(self$prefix)
      self$type    <- type
      self$step_name <- step_name
      self$keep_original_cols <- keep_original_cols

      if(length(self$columns) > 1 & self$type == "add") {
        rlang::abort("Add steps limit input columns to one.")
      }

      invisible(self)
    },
    # these are the base methods - can be overwritten in individual steps
    prep = function(new_data) {

      self$trained <- TRUE
      # self$result  <- vector(mode = "list", length(self$columns))
      #
      # names(self$result) <- file.path(self$id, self$columns, fsep = '_')

      invisible(self)
    },
    bake = function() {
      invisible(self)
    },
    tidy = function() {

      if (self$type == 'add') {

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
a <- function(..., x = 'blah', y = 1, z = NULL) {
  tmp <- c(`...` = enquos(...), environment())
  print(get_env(environment()))
  print(as.list(environment()))
}
#
# l <- function(..., x = 'blah', y = 1, z = NULL) {
#   b
# }
#
# z <- a(b, x = 'ada')
#
# do.call(l, z)



