#' #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' #
#' # Predict Regression Terms -----------------------------------------------------
#' #
#' #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' #' R6 Class
#' #'
#' #' `StepOlsResponse` Uses the Eigen C++ library fast versions to generate
#' #' predictions from different steps.
#' #'
#' #' @inheritParams Step
#' #'
#' #' @export
#' StepOlsResponse <- R6Class(
#'   classname = "step_ols_response",
#'   inherit = Step,
#'   public = list(
#'
#'     # step specific variables
#'     outcomes = NULL,
#'     predictors = NULL,
#'     coefficients = NULL,
#'     response_data = NULL,
#'     formula = NULL,
#'
#'     initialize = function(formula,
#'                           role = "augment",
#'                           ...) {
#'       # get function parameters to pass to parent
#'       # terms <- substitute(terms)
#'       env_list <- get_function_arguments()
#'       env_list$step_name <- "step_ols_response"
#'       env_list$type <- "supervise_augment"
#'       super$initialize(
#'         terms = NULL,
#'         env_list[names(env_list) != "terms"],
#'         ...
#'       )
#'
#'       self$formula <- formula
#'
#'       invisible(self)
#'     },
#'     bake = function(new_data, term_info, steps) {
#'
#'       x <- get_regression_data(new_data, term_info, id_type = "predictor")
#'       y <- get_regression_data(new_data, term_info, id_type = "outcome")
#'
#'
#'       self$coefficients <- determine_coefficients(x, y)
#'
#'       # column names in term info
#'       co_names <- x$term_info$variable
#'
#'       # print(x$term_info)
#'       resp <- list()
#'
#'       for (i in seq_along(steps)) {
#'         wh  <- collapse::whichv(x$term_info$step_index, i)
#'         co_name <- co_names[wh]
#'
#'         if (length(co_name) > 0) {
#'           co <- self$coefficients[wh, , drop = FALSE]
#'           resp[[i]] <- steps[[i]]$response(co)
#'           if (!"outcome" %in% names(resp[[i]])) {
#'             resp[[i]]$outcome <- rep(colnames(co), times = nrow(co))
#'           }
#'           if (!"term" %in% names(resp[[i]])) {
#'             resp[[i]]$term <- rep(co_name, times = ncol(co))
#'           }
#'         }
#'       }
#'
#'       resp <- collapse::rowbind(resp)
#'
#'       # save the response
#'       self$response_data <-
#'         append(resp, list(ols = rep.int(self$id, length(resp[[1]]))))
#'
#'       return(NULL)
#'     }
#'   )
#' )
