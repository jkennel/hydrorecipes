get_types <- function(data) {
  vapply(data, FUN = function(x) class(x)[1], FUN.VALUE = character(1))
}
