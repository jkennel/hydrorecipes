#' get_coefficients
#'
#' @param co coefficients from a model
#' @param step_info info on a recipe step
#'
#' @return the coefficients associated with a recipe step
#' @export
#'
get_coefficients <- function(co, step_info) {

  nms <- step_info$key
  co[which(names(co) %in% nms)]

}

#' response_lag
#'
#' @inheritParams get_coefficients
#'
#' @return a tibble containing the responses from lag terms
#' @export
#'
response_lag <- function(co, step_info) {

  co <- get_coefficients(co, step_info)

  tibble(co = co,
         cumulative = cumsum(co),
         lag = step_info$shift
  )

}

#' response_distributed_lag
#'
#' @inheritParams get_coefficients
#'
#' @return a tibble containing the responses from distributed lag terms
#' @export
#'
response_distributed_lag <- function(co, step_info) {

  knots <- step_info$knots
  bl <- basis_lag(knots)

  co <- get_coefficients(co, step_info)

  resp <- bl %*% co

  tibble(co = as.numeric(resp),
         cumulative = cumsum(as.numeric(resp)),
         lag = min(knots):max(knots))

}


#' response_harmonic
#'
#' @inheritParams get_coefficients
#'
#' @return a tibble containing the responses from harmonic terms
#' @export
#'
response_harmonic <- function(co, step_info) {

  co <- get_coefficients(co, step_info)
  sin_co <- co[grepl('sin',  names(co))]
  cos_co <- co[grepl('cos',  names(co))]

  tibble(frequency = step_info[grep('sin',  step_info$key),]$frequency,
         sin_co,
         cos_co,
         amplitude = sqrt(sin_co^2 + cos_co^2),
         phase = atan2(cos_co, sin_co))

}


get_step_and_coefficients <- function(rec) {
  rec$last_term_info
}


#' response
#'
#' @param fit a model fit object having a coefficients method.
#' @param rec a prepped recipe object..
#'
#' @return a list of response functions corresponding to each step
#' @export
#'
response <- function(fit, rec, data, ...) UseMethod("response")


#' @rdname predict_terms
#' @export
response.lm <- function(fit, rec) {

  co <- coefficients(fit)

  rec_steps <- tidy(rec)

  resp <- vector(mode = "list", length = nrow(rec_steps))
  names(resp) <- paste0(rec_steps$type, '_', rec_steps$step_name)

  for (i in 1:nrow(rec_steps)) {

    step_info <- tidy(rec, i)

    type <- rec_steps$type[i]

    if (type %in% c('lead', 'lag')) {
      out <- response_lag(co, step_info)
    } else if (type == 'distributed_lag') {
      out <- response_distributed_lag(co, step_info)
    } else if (type %in% c('harmonic', 'earthtide')) {
      out <- response_harmonic(co, step_info)
    } else {
      print(paste0('No response method for type: ', type))
      out <- NULL
    }

    if(!is.null(out)) {
      resp[[i]] <- out
    }

  }

  resp
}

#' @rdname predict_terms
#' @export
response.cv.glmnet <- function(fit, rec) {

  co <- coefficients(fit)
  co_names <- rownames(co)
  co <- as.vector(co)
  names(co) <- co_names

  rec_steps <- tidy(rec)

  resp <- vector(mode = "list", length = nrow(rec_steps))
  names(resp) <- paste0(rec_steps$type, '_', rec_steps$step_name)

  for (i in 1:nrow(rec_steps)) {

    step_info <- tidy(rec, i)

    type <- rec_steps$type[i]

    if (type %in% c('lead', 'lag')) {
      out <- response_lag(co, step_info)
    } else if (type == 'distributed_lag') {
      out <- response_distributed_lag(co, step_info)
    } else if (type %in% c('harmonic', 'earthtide')) {
      out <- response_harmonic(co, step_info)
    } else {
      print(paste0('No response method for type: ', type))
      out <- NULL
    }

    if(!is.null(out)) {
      resp[[i]] <- out
    }

  }

  resp
}

