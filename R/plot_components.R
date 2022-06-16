#' plot_components
#'
#' @param x the dataset
#' @param x_axis_column name of the column for the x coordinates
#' @param plot_columns names of the columns to plot
#' @param ... other arguments
#'
#' @return a ggplot2 plot
#' @export
#'
plot_components <- function(x,
                            x_axis_column = 'datetime',
                            plot_columns = NULL, ...) UseMethod("plot_components")

#' @export
plot_components.deconvolution <- function(x,
                                          x_axis_column,
                                          plot_columns, ...) {

  comp <- x[['components']]

  if(nrow(comp) > 0) {

      p <- plot_components(comp, x_axis_column, plot_columns)

  }
  p

}

#' @export
plot_components.data.frame <- function(x, x_axis_column = 'datetime', plot_columns = NULL, ...) {

  if(!is.null(plot_columns)) {
    plot_columns <- setdiff(plot_columns, x_axis_column)
  } else {
    plot_columns <- setdiff(names(x), x_axis_column)
  }

  nms <- unique(c(x_axis_column, plot_columns))
  x <- x[, nms]
  x <- x |> pivot_longer(!starts_with(x_axis_column))

  x$name <- factor(x$name, levels = plot_columns, labels = plot_columns)

  p <- ggplot(x, aes_string(x = x_axis_column, y = 'value')) +
    geom_line() +
    facet_wrap(name~., scales = 'free_y', nrow = length(unique(x[['name']]))) +
    theme_bw()

  p
}
