#' pad_num
#'
#' This function creates a sequence of numbers with 0 padding based on the series
#' length.
#'
#' @param n the number of columns
#' @param pad the character to use for the padding.
#'
#' @return
#' @export
#'
pad_num <- function(n, pad = "0") {

  width <- floor(log10(n)) + 1L

  formatC(seq_len(n),
          width = width,
          format = "d",
          flag = "0")

}

name_columns <- function(id, column_name, n) {
  if (is.null(column_name)) {
    if (n < 2) {
      return(file.path(id, fsep = '_'))
    }
    return(file.path(id, pad_num(n), fsep = '_'))
  }
  if (n < 2) {
    return(file.path(id, column_name, fsep = '_'))
  }

  file.path(id, column_name, pad_num(n), fsep = '_')

}




#' Make a random identification field for steps
#'
#' @export
#' @param prefix A single character string
#' @param len An integer for the number of random characters
#' @return A character string with the prefix and random letters separated by
#'  and underscore.
#'
#' @useDynLib frecipes, .registration = TRUE
#' @importFrom R6 R6Class
#' @importFrom Rcpp sourceCpp
#' @keywords internal
rand_id <- function(prefix = "step", len = 5L) {
  candidates <- c(letters, LETTERS, paste(0:9))
  paste(prefix,
        paste0(sample(candidates, len, replace = TRUE), collapse = ""),
        sep = "_"
  )
}

# parse_variables <- function(env, arg_nms) {
#
#       env_nms <- names(env)
#       super_nms <- arg_nms
#       env_nms_sub <- intersect(super_nms, env_nms)[-1L]
#
#       inputs <- c(
#         as.list(substitute(test, environment())),
#         as.list(environment())[env_nms_sub]
#       )
#
# }

