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
    new_columns = NULL,
    step_name = NULL,
    keep_original_cols = TRUE,
    id = NULL,
    prefix = NULL,
    result = NULL,


    varying = NULL, # list(name = , initial = , lower = , upper = )
    rerun = TRUE,
    n_na_max = NULL,

    check = NULL,



    initialize = function(terms, ...) {

      if (!missing(terms)) {
        if (length(terms) == 1L) {
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

      # don't want flattened list so use list(...)
      self$varying <- list(...)[["varying"]]

      # super specific values
      if (is.null(self$prefix)) {
        self$prefix <- gsub("step_", "", self$step_name)
      }
      self$id <- rand_id(self$prefix)

      invisible(self)
    },

    # these are the base methods - can be overwritten in individual steps
    prep = function(column_names) {
      # self$columns <- get_terms_from_info(self$terms, column_names)
      self$trained <- TRUE

      invisible(self)
    },
    bake = function(s) {
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
        x = rep.int(NA_real_, n),
        variable = rep.int("coefficient", n),
        value = as.vector(co),
        step_id = rep.int(self$id, n),
        outcome = rep(colnames(co), each = n_each),
        term = rep.int(self$prefix, n),
        step_columns = rep.int(paste(self$columns, collapse = "_"), n)
      )

    },
    get_fields = function() {
      sapply(self, class)
    },
    get_number_columns = function() {
      length(self$new_columns)
    },
    get_result = function(column_name = NULL) {

      if (is.null(column_name)) {
        return(self$result)
      }

      nms <- names(self$result)
      if (column_name %iin% nms) {
        return(self$result[column_name])
      }

      return(NULL)
    },

    is_coef = function() {

      if (!is.null(self[["varying"]])) {
        if(all(self[["varying"]][["name"]] == "coef")) {
          return(TRUE)
        }
      }

      return(FALSE)
    },

    set_result = function(values) {
      self$result <- values

      return(self)
    },
    update_varying = function() {

      if (!is.null(self[["varying"]])) {
        if (self$is_coef()) {
          n <- self$get_number_columns()
          if (length(self[["varying"]][["name"]]) == 1L & n > 1L) {
            self[["varying"]] <- list(name  = rep.int("coef", n),
                                      start = rep.int(self[["varying"]][["start"]], n),
                                      lower = rep.int(self[["varying"]][["lower"]], n),
                                      upper = rep.int(self[["varying"]][["upper"]], n))
          }
        }
      }

    },
    update_step = function(field_name, field_value) {

      for (i in seq_along(field_name)) {

        if(field_name[[i]] %iin% names(self)) {
          self[[field_name[i]]] <- field_value[i]
        }

      }

      self$rerun <- TRUE

    }

  )
)




