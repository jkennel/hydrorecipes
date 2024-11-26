#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Building Block for a Recipe --------------------------------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Step <- R6Class(
  classname = "step",
  public = list(
    type = NULL, # check, add, remove, update/modify

    # base steps
    terms = NULL,
    role = NULL,
    trained = FALSE,
    skip = FALSE,
    columns = NULL,
    step_name = NULL,
    keep_original_cols = TRUE,
    id = NULL,
    prefix = NULL,
    result = NULL,


    varying = NULL, # list(name = , initial = , lower = , upper = )
    rerun = TRUE,

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

      self$varying <- dots$varying

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
      n_each = nrow(co)
      n <- length(co)

      list(
        x = rep(NA_real_, n),
        variable = rep("coefficient", n),
        value = as.vector(co),
        step_id = rep(self$id, n),
        outcome = rep(colnames(co), each = n_each)
      )

    },

    get_fields = function() {
      sapply(self, class)
    },
    get_result = function(column_name = NULL) {

      if (is.null(column_name)) {
        return(self$result)
      }

      nms <- names(self$result)
      if (column_name %in% nms) {
        return(self$result[column_name])
      }

      return(NULL)
    },
    set_result = function(values) {
      self$result <- values

      return(self)
    },
    update_step = function(field_name, field_value) {

      self[[field_name]] <- field_value
      self$rerun <- TRUE

    }

  )
)




