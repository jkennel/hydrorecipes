create_new_column_names <- function(self, single = TRUE) {

  if(single) {

    self$new_columns <- file.path(
      self$prefix,
      self$columns
    )

  } else {

    self$new_columns <-
      file.path(
        self$prefix,
        rep(self$columns, length(self$names_vector)),
        "_",
        self$names_vector
      )
  }

  invisible(self)
}
