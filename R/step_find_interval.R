#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Divide a Term into Intervals and do Dummy Encoding ---------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' R6 Class
#'
#' `StepFindInterval` divides a series into intervals and then
#' performs dummy encoding.
#'
#' @param vec a vector of break points
#'
#' @inheritParams Step
#'
#' @export
StepFindInterval <- R6Class(
  classname = "step_find_interval",
  inherit = Step,
  public = list(

    # step specific variables
    vec = NULL,
    n_vec = NULL,

    #' @description
    #' @inheritParams StepAddVars
    #' @return A new `Step`.
    initialize = function(terms,
                          vec,
                          role = "predictor",
                          ...) {
      # get function parameters to pass to parent
      terms <- substitute(terms)
      env_list <- get_function_arguments()
      env_list$step_name <- "step_find_interval"
      env_list$type <- "add"
      super$initialize(
        terms = terms,
        env_list[names(env_list) != "terms"],
        ...
      )



      # step specific values
      self$vec <- sort(vec)
      self$n_vec <- length(vec)

      invisible(self)
    },
    bake = function(new_data) {
      column_name <- self$columns

      self$new_columns <- c()

      dum <- list()
      for (i in seq_along(column_name)) {
        dum[[i]] <- to_dummy_list(unclass(new_data)[[i]], self$vec)

        nn <- name_columns(
          self$prefix,
          column_name[i],
          length(dum[[i]])
        )

        names(dum[[i]]) <- nn
        self$new_columns <- c(self$new_columns, nn)
      }
      unlist(dum, recursive = FALSE)
    }
  )
)

# library(R6)
# library(recipes)
# library(hydrorecipes)
# set.seed(123)
# x <- sort(rnorm(5e7))
# y <- x
# nd <- list(x=x, y = y)
# bench::mark(
#   a <- StepFindInterval$new(x, vec = -3:3)$prep()$bake(new_data = nd),
#   b <- StepFindInterval$new(x, vec = -3:3)$prep()$bake2(new_data = nd),
#   # b_spline(x, -7:7),
#   iterations = 10,
#   check = FALSE
# )
# all.equal(a,b, check.attributes = FALSE)
#
# bench::mark(
#   tmp <- recipe(y~x, data = tibble(x=x, y = y)) |>
#     step_find_interval(x, y, vec = -7:7) |>
#     prep() |>
#     bake(new_data = NULL),
#   iterations = 5
# )
#
#
#
#
# n <- 1000000
# n_cols <- 3
# v <- sample(1:n_cols, size = n, replace = TRUE)
# f <- factor(v)
# y <- data.frame(z = 1:n, x = f,x2 = f, v = v, v2 = v)
#
# m <- tt2$steps[[1]]$result[[1]]
#
#
# # dummy_to_factor <- factor(m %*% 1:ncol(m), labels = colnames(m))
#
# microbenchmark::microbenchmark(
#   tt0 <- Recipe$new(z~x, data = y)$
#     add_step(StepFindInterval$new(x, vec = c(2,8)))$
#     add_step(StepFindInterval$new(x, vec = c(5,7)))$
#     prep()$
#     bake()$
#     bc(),
#   tt0 <- Recipe$new(z~x, data = y)$
#     add_step(StepFindInterval$new(x, vec = c(2,8), one_hot = FALSE))$
#     add_step(StepFindInterval$new(x, vec = c(5,7), one_hot = FALSE))$
#     prep()$
#     bake(),
#
#   # tt1 <- Recipe$new(z~x, data = y)$
#   #   step_find_interval(x, vec = c(2,8))$
#   #   step_find_interval(x, vec = c(5,7))$
#   #   prep()$
#   #   bake()$
#   #   bc(),
#   #
#   # # tt2 <- Recipe$new(z~x, data = y)$
#   # #   step_find_interval(x, vec = c(2,8), one_hot = FALSE)$
#   # #   prep()$
#   # #   bake(y),
#   #
#   tt2_5 <- frecipe(z~x, data = y) |>
#     fstep_find_interval(x, vec = c(2,8)) |>
#     fprep() |>
#     fbake(y),
#   # tt2_6 <- frecipe(z~x, data = y) |>
#   #   fstep_find_interval(x, vec = c(2,8), one_hot = TRUE) |>
#   #   fprep() |>
#   #   fbake(y),
#
#   tt3 <- recipe(z~x+x2, data = y) |>
#     step_dummy(x)|>
#     step_dummy(x2)|>
#     prep() |>
#     bake(new_data = y),
#
#   tt4 <- recipe(z~v+v2, data = y) |>
#     step_cut(v, breaks = c(2,8))|>
#     step_cut(v2, breaks = c(2,8))|>
#     prep() |>
#     bake(new_data = NULL),
#
#   times = 2
#
# )
#
# system.time({
#   tt <- StepFindInterval$new(x, vec = 1:10)
#   tt$prep()
#   tt$bake(new_data = y)
# })
#
# system.time({
#   tt <- StepFindInterval$new(x, vec = 1:10, one_hot = FALSE)
#   tt$prep()
#   tt$bake(new_data = y)
# })
#
# system.time({
#   r <- recipe(~x, data = y) |>
#     step_dummy(x) |>
#     prep() |>
#     bake(new_data = NULL)
# })
#
# system.time({
#   r <- recipe(~v, data = y) |>
#     step_cut(v, breaks = 1:10) |>
#     prep() |>
#     bake(new_data = NULL)
# })
#
#
# tt$bake()
# n <- 1e1
# x <- data.table(x = 1:n)
# microbenchmark::microbenchmark(
#   min(unclass(x)[[1]]),
#   max(unclass(x)[[1]]),
#   range(unclass(x)[[1]])
# )
#
#
# tmp <- function(...) {
#   (rlang::as_name(enquo(...)))
# }
#
# fo <- function(...) {
#   d <- list(...)
#   print(rlang::quo_text(enquo(d)))
#   print(rlang::as_name(enquo(d)))
#   print(rlang::quo_name(enquo(d)))
# }
#
# a <- fo(b)
#
#
# tmp <- function(..., df) {
#   vapply(ensyms(...),
#          FUN = rlang::as_name,
#          FUN.VALUE = character(1))}
#
#
#
# dt <- data.table(a = 1:10, b = 1:10, df = dt)
# bench::mark(
#
#   d <- tmp(e, b, df = dt),
#   d2 <- tmp2(e, b, df = dt)
#
# )
#
# bench::mark(
#   paste0('a', 'b', collapse = '', recycle0 = FALSE),
#   file.path('a', 'b', fsep = '')
# )
#
# n <- 1000000
# a <- (data.frame(x = 1:n))
# b <- matrix(c(1:n, 1:n), ncol = 2)
# colnames(b) <- c('a', 'b')
# d <- list(a =1:n, b = 1:n)
# bench::mark(
#   q <- bind_cols(a, b,b, .name_repair = 'minimal'),
#   r <- bind_cols(a,d, d, .name_repair = 'minimal'),
# )
