#' plot_response
#'
#' @param x the dataset
#' @param titles the plot titles - one for each response
#' @param ... other arguments
#'
#' @return list of ggplot2 plots
#' @export
#'
plot_response <- function(x, titles = NULL, ...) UseMethod("plot_response")

#' @export
plot_response.deconvolution <- function(x, titles = NULL, ...) {

  # remove entries with no response
  resp <- x[['response']]
  resp <- resp[!sapply(resp, is.null)]

  # check titles
  if (!is.null(titles)) {
    if(length(titles) != length(resp)) {
      stop('title length and number of responses should be equal')
    }
  }

  # generate plots
  if(length(resp) > 0) {

    p_list <- list()
    for(i in seq_along(resp)) {
      if(is.null(titles)) {
        title <- names(resp)[i]
      } else {
        title <- titles[i]
      }
      p_list[[i]] <- plot_response(resp[[i]], titles = title)
    }

  }
  p_list
}


#' @export
plot_response.data.frame <- function(x, titles = NULL, ...) {
  frequency <- amplitude <- phase <- value <- lag <- cumulative <- NULL
  nms <- names(x)
  if('lag' %in% nms) {
    p <- ggplot(x, aes(x = lag, y = cumulative)) +
      geom_line() +
      ggtitle(titles) +
      scale_y_continuous(limits = c(0, max(x$cumulative, na.rm = TRUE))) +
      theme_bw()
  } else if('frequency' %in% nms) {
    x <- x[, c('frequency', 'amplitude', 'phase')]
    x <- x |> pivot_longer(!frequency)
    p <- ggplot(x, aes(x = frequency, xend = frequency, y = 0, yend = value)) +
      geom_segment() +
      ggtitle(titles) +
      facet_wrap(name~., scales = 'free_y', nrow = 2) +
      theme_bw()
  }
  p
}
