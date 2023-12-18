#' R6 Class
#'
#' `Recipe` holds the steps to create model features
#'
#' @importFrom recipes recipes_eval_select
#'
#' @export
Recipe <- R6Class(

  classname = 'recipe',

  public = list(

    #' @field the model formula.
    formula   = NULL,
    #' @field information about original variables.
    var_info  = NULL,
    #' @field information about original and predicted variable.
    term_info = NULL,
    #' @field information about original and predicted variable.
    last_term_info = NULL,
    #' @field the recipe steps.
    steps     = list(),
    #' @field the data.
    template  = list(),
    #' @field variable factor levels.
    levels    = NULL,
    #' @field variables that are retained.
    retained  = NA,
    #' @field packages required
    requirements = NULL,
    #' @field training info
    tr_info = NULL,
    #' @field original factor levels
    orig_lvls = NULL,


    #' @field whether the model is trained.
    trained = NULL,

    #' @field list that holds the created model features.
    result = list(),

    #' @field list that holds the checks.
    checks = list(),

    vars = NULL,

    #' @description
    #' Create a new recipe object.
    #' @param formula The model formula. It cannot contain operations.
    #' @param data list, data.frame, data.table, tibble of data. They will all
    #' be treated as lists.
    #' @return A new `Recipe` object.
    initialize = function(formula, data) {

      # specify data used with formula notation
      self$formula <- formula

      # parse the formula
      vars_list <- get_formula_vars_2(formula = formula, data = data)
      self$vars <- unlist(vars_list, use.names = FALSE)
      self$template <- unclass(data)[unique(self$vars)]

      # variable info
      self$var_info  <- list(
        variable = self$vars,
        type = get_types(self$template),
        roles = rep.int(x = c('predictor', 'outcome'),
                        times = vapply(vars_list, FUN = length, FUN.VALUE = numeric(1))),
        source = rep.int('original', times = length(self$vars)),
        step_index = rep.int(0L, times = length(self$vars)),
        step_name = rep.int("initial", times = length(self$vars)))

      self$term_info <- self$var_info


      self$requirements <- list(
        bake = setNames(object = logical(), nm = character())
      )

      # self$add_step(StepAddVars$new(vars = vars))#$prep()$bake(data)

      invisible(self)

    },

    #' @description
    #' Add a step to the recipe.
    #' @param step The step to add.
    #' @return An updated `Recipe` object.
    add_step = function(step) {

      step$columns = recipes::recipes_eval_select(step$terms,
                                                  self$template,
                                                  qDF(self$var_info))


      self$steps <- append(self$steps, step)
      invisible(self)

    },

    #' @description
    #' Do prep operations.
    #' @return An updated `Recipe` object.
    prep = function(retain = TRUE) {
      self$tr_info <- self$train_info(self$template)

      for (i in seq_along(self$steps)) {
        columns <- self$steps[[i]]$columns
        self$steps[[i]]$prep(unclass(self$template)[columns])
      }

      self$retained <- retain
      invisible(self)
    },

    #' @description
    #' Create the dataset.
    #' @param new_data The input data to the recipe.
    #' @return An updated `Recipe` object with a result that holds a list of
    #' features.
    bake = function(new_data = NULL, type = "list") {
      types <- self$get_step_types()

      if(is.null(new_data)) {
        self$result <- self$template
      } else {
        self$result <- unclass(new_data)[unique(self$vars)]
      }

      for (i in seq_along(types)) {
        columns <- self$steps[[i]]$columns

        # only one column allowed
        if (types[i] == "add"){

          self$result <- append(self$result,
                                self$steps[[i]]$bake(unclass(self$result)[[columns]]))

        } else if (types[i] == "modify") {
          self$result <- modifyList(self$result,
                                    self$steps[[i]]$bake(unclass(self$result)[columns]))

        } else if (types[i] == "supervised_add") {
          self$result <- append(self$result,
                                self$steps[[i]]$bake(unclass(self$result), self$term_info))

        } else if (types[i] == "check") {
          self$checks <- append(self$checks,
                                self$steps[[i]]$bake(unclass(self$result)[[columns]]))
        }
        self$update_term_info(step_name = self$steps[[i]]$step_name, step_index=i)
      }


      invisible(self)

    },
    ## short summary of training set.
    train_info = function(x) {
      data.frame(
        nrows = length(x[[1L]]),
        ncomplete = collapse::fsum(!collapse::missing_cases(x))
      )
    },

    update_term_info = function(source = "derived",
                                roles = "predictor",
                                type = "numeric",
                                step_name,
                                step_index) {

      variable <- setdiff(names(self$result), self$term_info$variable)

      n <- length(variable)

      self$term_info$variable <- c(self$term_info$variable, variable)
      self$term_info$roles <- c(self$term_info$roles,
                                rep.int("predictor", times = n))
      self$term_info$source <- c(self$term_info$source,
                                 rep.int(source, times = n))
      self$term_info$type <- c(self$term_info$type,
                               rep.int(type, times = n))
      self$term_info$step_index <- c(self$term_info$step_index,
                                     rep.int(step_index, times = n))
      self$term_info$step_name <- c(self$term_info$step_name,
                                    rep.int(step_name, times = n))
    },

    get_step_types = function() {
      vapply(self$steps, FUN = function(x) x$type, FUN.VALUE = character(1))
    },

    data = function(type = "df") {

      if (type == "df") {
        return(collapse::qDF(self$result))
      }

      if (type == "m") {
        return(collapse::qM(self$result))
      }

      if(type == "dt") {
        return(collapse::qDT(self$result))
      }

      if(type == "tbl") {
        return(collapse::qTBL(self$result))
      }

      return(self$result)
    }


  )
)




# library(tibble)
# library(data.table)
# x <- rnorm(1e6)
# xl <- list(x)
# y <- x
# df <- as.data.frame(x)
# dt <- as.data.table(x)
# tb <- as_tibble(x)
# bench::mark(
#   dt[, y := y],
#   dt$y <- y,
#   df$y <- y,
#   df[["y"]] <- y,
#   add_column(tb, y, .name_repair = "minimal"),
#   set(dt, j = "y", value = y),
#   xl[["y"]] <- y,
#   xl[[1]] <- y,
#   check = FALSE
# )
#
# bench::mark(
#   as.list(df)[names(df)],
#   unlist(df, recursive = FALSE)[names(df)],
#   df[names(df)],
#   check = FALSE
# )
#
# n <- 10
# a = 1:n
# b = rnorm(n)
# bench::mark(
#   data.frame(a, b),
#   tibble(a,b),
#   data.table(a,b),
#   list(a, b),
#   qDF(list(a, b)),
#   check = FALSE
# )
#
#
# fstep_find_interval <- function(rec,
#                                 ...,
#                                 vec,
#                                 one_hot = TRUE,
#                                 role = "predictor",
#                                 skip = FALSE,
#                                 prefix = "step_find_interval_",
#                                 id = rand_id("step_find_interval"),
#                                 keep_original_cols = FALSE) {
#
#   rec$add_step(step = StepFindInterval$new(...,
#                                            vec = vec,
#                                            one_hot = one_hot,
#                                            role = role,
#                                            skip = skip,
#                                            prefix = prefix,
#                                            id = id,
#                                            keep_original_cols = keep_original_cols))
#
#   invisible(rec)
#
# }
#
# # compatible with recipe format
# frecipe <- function(formula, data) {
#   Recipe$new(formula, data)
# }
#
# fprep <- function(rec) {
#   rec$prep()
# }
#
# fbake <- function(rec, new_data = NULL) {
#   rec$bake(new_data)
# }
#
# m <- matrix(1:10000000, ncol = 10)
# colnames(m, paste0("n", 1:10))
# tmp <- list(a = 1:1000000, b = 1:1000000, d = 1:1000000, e = 1:1000000, df = data.frame(z = 1:1e6, y = rnorm(1e6)))
#
# tt <- qDF(tmp)
#
# bench::mark(qM(tmp),
#             bind_cols(tmp),
#             qDF(tmp),
#             qM(qDF(tmp)),
#             do.call(cbind, tmp),
#             # tt['b'] <- matrix(1:90000, ncol = 9),
#             cbind2(tmp), check = FALSE)
#
# dat <- matrix(1:2e6, ncol = 2)
# microbenchmark::microbenchmark(
#   as_tibble(dat),
#   as.data.frame(dat),
#   as.data.table(dat),
#   qDT(dat)
# )
#
#
# library(recipes)
# library(data.table)
# library(R6)
# library(microbenchmark)
# library(collapse)
#
# n <- 1000000
# dat <- data.frame(x = 1:n,
#                   y = 1:n,
#                   z = 1:n,
#                   q = rnorm(n),
#                   r = qF(rnorm(n)),
#                   a = qF(sample(1:10, n, replace = TRUE)),
#                   b = qF(sample(1:3, n, replace = TRUE)))
#
# formula <- as.formula('r~b+ a + q')
# bench::mark(
#   (recipe(formula = formula, data = dat)),
#   # (recipe(x = dat, vars = c('y', 'x'))),
#   tmp <- (Recipe$new(formula = formula, data = dat)),
#   check = FALSE
# )
# # microbenchmark::microbenchmark(
# #   frec <- Recipe$new(formula = formula, data = dat),
# #   rec <- recipe(formula = formula, data = dat),
# #   recipe(x = dat,
# #          vars = c('y', 'x')),
# #   times = 10
# # )
#
# # m <- as.matrix(dat)
#
# a <- rlang::f_rhs(formula)
# vapply(a,
#        FUN = function(x) paste0(deparse(x)),
#        FUN.VALUE = character(1))[-1L]
#
# microbenchmark::microbenchmark(
#   attr(model.frame(formula, dat, subset = 1L), "terms"),
#   model.frame(formula, dat, subset = 1L),
#   # terms(formula, dat[1,]),
#   attr(model.frame(formula, dat[1,]), "terms"),
#   rlang::f_lhs(formula)
# )
#
# formula <- as.formula('.~.')
#
#
# get_rhs_vars_2(formula, dat[1,])
# get_rhs_vars(formula, dat[1,])
#
#
# microbenchmark(
#   get_formula_vars(formula, dat),
#   get_formula_vars_2(formula, dat),
#   # get_rhs_vars_2(formula, dat, no_lhs = FALSE),
#   # get_lhs_vars_2(formula, dat),
#   parse_formula(formula),
#   parse_formula_2(formula)
#   # form2args(formula, dat)
#
# )
#
# bench::mark(colnames(dat), names(dat))
# bench::mark(
#   get_lhs_vars(formula, data),
#   get_lhs_vars_2(formula, dat),
# )
#
#
#
# require(microbenchmark)
# require(data.table)
# k <- "keycol"
# N <- 1e7
# DT = data.table(a = runif(N), b = rnorm(N))
# DF = data.frame(a = runif(N), b = rnorm(N))
# L <- list(a = runif(N), b = rnorm(N))
# sl <- seq_len(nrow(DT))
# tb <- tibble(a = runif(N), b = rnorm(N))
# ans <- capture.output(microbenchmark(
#     DT[,keycol := sl],
#     DT$keycol <- sl,     #as mentioned in vignette, this is slow
#     DT[["keycol"]] <- sl,
#     DT[,"keycol"] <- sl,
#     L[[k]] <- sl,
#     L[["keycol"]] <- sl,
#     L$keycol <- sl,
#     DF$keycol <- sl,
#     DF[["keycol"]] <- sl,
#     DF[,"keycol"] <- sl,
#
#     times = 10L))
# message(paste0("#",ans,"\n"))
#
# vars <- c("a", "b", "b")
# ans <- capture.output(bench::mark(
#     L[vars],
#     tb[, vars],
#     check = FALSE))
# message(paste0("#",ans,"\n"))
