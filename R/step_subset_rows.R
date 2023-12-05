#' R6 Class
#'
#' `StepSubsetRows` selects rows from output.
#' @inheritParams Step
#'
#' @export
StepSubsetRows <- R6Class(
  classname = 'step_subset_rows',
  inherit = Step,

  public = list(
    row_numbers = NULL,
    # step specific variables
    initialize = function(...,
                          row_numbers,
                          role = "modify",
                          skip = FALSE,
                          keep_original_cols = FALSE) {

      # get function parameters to pass to parent
      step_name    <- "step_subset_rows"
      type         <- 'modify'
      inputs <- c(
        as.list(rlang::quos(...)),
        rlang::env_get_list(env = environment(),
                            formalArgs(super$initialize)[-1L])
      )
      do.call(super$initialize, inputs)

      self$row_numbers <- as.integer(row_numbers)
      invisible(self)
    },
    bake = function(new_data) {
      unclass(collapse::qDT(new_data)[self$row_numbers,])
    }
  )
)

# library(data.table)
# n <- 1e7
# df <- data.frame(age=sample(1:65,1e7,replace=TRUE),x=rnorm(1e7),y=rpois(1e7,25))
# dt <- as.data.table(df)
# vec <- df[,1]
#
# subsetter <- function(x, y) {
#   lapply(unlist(x, recursive = FALSE), '[', y)
# }
# subsetter3 <- function(x, y) {
#   lapply(unlist(x, recursive = FALSE), vctrs::vec_slice, y)
# }
# subsetter2 <- function(x, y) {
#   lapply(x, '[', y)
# }
# dt <- qDT(unlist(tmp$result, recursive = FALSE))
# df <- qDF(unlist(tmp$result, recursive = FALSE))
# tbl <- qTBL(unlist(tmp$result, recursive = FALSE))
# m <- qM(unlist(tmp$result, recursive = FALSE))
#
# to_rem <- as.integer(seq(1L, 1e7L, 2L))
# w <- unlist(tmp$result, recursive = FALSE)
# class(w) <- "data.table"
# bench::mark(
# #   aa<-subsetter(tmp$result, to_rem),
# #   aa<-subsetter3(tmp$result, to_rem),
# #   bb<-subsetter2(df, to_rem),
#   ww <- unclass(qDT(unlist(tmp$result, recursive = FALSE))[to_rem]),
#   m[to_rem,],
#   w[to_rem],
#   dt[to_rem],
#   dt[to_rem,],
#   tbl[to_rem,],
#   check = FALSE,
#   min_iterations = 10
# )
