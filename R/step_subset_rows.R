#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Subset dataset rows ----------------------------------------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
StepSubsetRows <- R6Class(
  classname = "step_subset_rows",
  inherit = Step,
  public = list(
    # step specific variables
    row_numbers = NULL,
    initialize = function(row_numbers = 1L,
                          role = "modify",
                          ...) {
      # get function parameters to pass to parent
      env_list <- get_function_arguments()
      env_list$step_name <- "step_subset_rows"
      env_list$type <- "model"

      super$initialize(
        terms = NULL,
        env_list[names(env_list) != "terms"],
        ...
      )

      self$row_numbers <- as.integer(row_numbers)
      invisible(self)

    },
    bake = function(r) {

      self$result <- collapse::ss(r$get_result(type = "list"),
                                  i = self$row_numbers, check = FALSE)

      r$template_step <- length(r$time_bake) + 1L

      return(NULL)
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
