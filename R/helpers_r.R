pad_num <- function(n, pad = "0") {

  width <- floor(log10(n)) + 1L

  formatC(seq_len(n),
          width = width,
          format = "d",
          flag = "0")

}



create_new_column_names <- function(self, len = 1) {

  new_columns <-
    file.path(
      self$id,
      pad_num(len),
      fsep = "_"
    )

  new_columns
}

