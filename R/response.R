#' get_coefficients
#'
#' @description Get the coefficients by name that are related to a particular step_*
#'
#' @param co coefficients from a model
#' @param step_info info on a recipe step
#'
#' @return the coefficients associated with a recipe step
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
#'
response_lag <- function(co, step_info) {

  co <- get_coefficients(co, step_info)

  tibble(coefficient = co,
         cumulative = cumsum(co),
         x = step_info$shift
  )

}


#' response_distributed_lag
#'
#' @inheritParams get_coefficients
#'
#' @return a tibble containing the responses from distributed lag terms
#'
response_distributed_lag <- function(co, step_info) {

  knots <- step_info$knots
  bl <- basis_lag(knots)

  co <- get_coefficients(co, step_info)

  resp <- bl %*% co

  tibble(coefficient = as.numeric(resp),
         cumulative = cumsum(as.numeric(resp)),
         x = min(knots):max(knots))

}


#' response_harmonic
#'
#' @inheritParams get_coefficients
#'
#' @return a tibble containing the responses from harmonic terms
#'
response_harmonic <- function(co, step_info) {

  co <- get_coefficients(co, step_info)
  sin_coefficient <- co[grepl('sin',  names(co))]
  cos_coefficient <- co[grepl('cos',  names(co))]

  tibble(x = step_info[grep('sin',  step_info$key),]$frequency,
         sin_coefficient,
         cos_coefficient,
         amplitude = sqrt(sin_coefficient^2 + cos_coefficient^2),
         phase = atan2(cos_coefficient, sin_coefficient))

}


#' response
#'
#' @description This function takes a model object and extracts the responses
#'  from step_distributed_lag, step_lag, step_lead, step_lead_lag, step_harmonic
#'  and step_earthtide
#'
#' @param fit a model fit object having a coefficients method.
#' @param rec a prepped recipe object.
#' @param verbose print names of steps with no response methods
#' @param ... currently not used
#'
#' @return a list of response functions corresponding to each step
#' @export
#'
response <- function(fit, rec, verbose = FALSE, ...) UseMethod("response")


#' @rdname response
#' @export
response.lm <- function(fit, rec, verbose = FALSE, ...) {

  co <- coefficients(fit)

  response.numeric(co, rec, verbose = FALSE,...)

}


#' @rdname response
#' @export
response.cv.glmnet <- function(fit, rec, verbose = FALSE, ...) {

  co <- coefficients(fit)
  co_names <- rownames(co)
  co <- as.vector(co)
  names(co) <- co_names

  response.numeric(co, rec, ...)
}

#' @rdname response
#' @export
response.numeric <- function(fit, rec, verbose = FALSE, ...) {

  x <- NULL

  rec_steps <- tidy2(rec)

  resp <- vector(mode = "list", length = nrow(rec_steps))
  names(resp) <- paste0(rec_steps$type, '_', rec_steps$step_name)

  for (i in 1:nrow(rec_steps)) {

    step_info <- tidy2(rec, i)

    type <- rec_steps$type[i]

    if (type %in% c('lead', 'lag', 'lead_lag')) {
      out <- response_lag(fit, step_info)
    } else if (type == 'distributed_lag') {
      out <- response_distributed_lag(fit, step_info)
    } else if (type %in% c('harmonic', 'earthtide')) {
      out <- response_harmonic(fit, step_info)
    } else {
      if(verbose) {
        print(paste0('No response method for type: ', type))
      }
      out <- NULL
    }


    if(!is.null(out)) {

      out <- pivot_longer(out, cols = !x)
      out$type <- type
      out$term <- unique(step_info$terms)

      resp[[i]] <- out

    }

  }

  bind_rows(resp)
}
