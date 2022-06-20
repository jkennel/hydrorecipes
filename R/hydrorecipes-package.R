#' hydrorecipes: Hydrogeology Steps for the 'recipes' Package
#'
#' The package provides a few steps for the analysis of water level
#' datasets of moderately size (millions of rows) using the `recipes` package.
#' It aims to work seamlessly within the `tidymodels` framework of packages.
#' While the package was designed for water levels and pore-pressures, the
#' steps may be useful in other fields of study. Examples in the introduction
#' focus on barometric pressure and Earth tides, however, the steps may also be
#' applied to precipitation, river stage, ocean tides, and pumping responses.
#'
#'
#'
#' You can learn about the hydrorecipes package in the vignettes:
#' \code{browseVignettes(package = "hydrorecipes")}
#' @keywords internal
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
#' @useDynLib hydrorecipes, .registration = TRUE
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL
