# #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# #
# # Add variables to Recipe Step -------------------------------------------------
# #
# #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# StepDataTransducer <- R6Class(
#
#   classname = "step_data_transducer",
#   inherit = Step,
#
#   public = list(
#
#     file_name = NA_character_,
#     file_type = NA_character_,
#     dots = NULL,
#     data = NULL,
#
#     # step specific variables
#     initialize = function(file_name,
#                           file_type = "RBR",
#                           role = "data",
#                           ...) {
#
#       # get function parameters to pass to parent - strings can be passed
#       terms <- substitute(terms)
#       env_list <- get_function_arguments()
#       env_list$step_name <- "step_data_transducer"
#       env_list$type <- "data"
#       super$initialize(
#         terms = terms,
#         env_list[names(env_list) != "terms"],
#         ...
#       )
#       self$dots <- list(...)
#       self$file_name <- file_name
#       self$file_type <- file_type
#
#
#       invisible(self)
#     },
#
#     bake = function() {
#
#       dat <- list()
#       for (i in seq_along(self$file_name)) {
#         if (self$file_type == "RBR") {
#           self$data[[i]] = rsk::Rsk$new(self$file_name[[i]], ...)
#         }
#         if (self$file_type == "diver") {
#           # self$data[[i]] = rsk::Rsk$new(self$file_name[[i]], ...)
#         }
#         if (self$file_type == "levellogger") {
#           # self$data[[i]] = rsk::Rsk$new(self$file_name[[i]], ...)
#         }
#
#
#         # serial <- self$data[[i]]$serial
#         # model <- self$data[[i]]$instruments$model
#         d <- melt(self$data[[i]], id.vars = "datetime")
#         d$file_name <- basename(self$file_name[i])
#         dat[[i]] <- d
#       }
#
#
#       # combine
#       dat <- rbindlist(dat)
#
#       return(unclass(self$data))
#
#     }
#
#   )
# )
#
