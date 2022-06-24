#' get_coefficients
#'
#' @description Get the coefficients by name that are related to a particular step_*
#'
#' @param co named coefficients from a model fit
#' @param step_info info on a `recipe` step
#'
#' @return the coefficients associated with a particular `recipe` `step`
#'
#' @noRd
get_coefficients <- function(co, step_info) {

  nms <- step_info$key

  co_sub  <- co[which(names(co) %in% nms)]

  # handle cases where numbers are not present in the key
  if(length(co_sub) == 0) {
    co_nms <- names(co)
    co_nms <- sub(nms, '', co_nms)
    co_nms <- sub('_', '', co_nms)
    co_sub <- co[!suppressWarnings(is.na(as.numeric(co_nms)))]
  }
  co_sub[is.na(co_sub)] <- 0

  co_sub

}


#' response_lag
#'
#' @description Calculate the response for a lagged set
#'
#' @inheritParams get_coefficients
#'
#' @return a `tibble` containing the responses from lag terms
#'
#' @noRd
#'
response_lag <- function(co, step_info) {

  co <- get_coefficients(co, step_info)

  tibble(coefficient = co,
         cumulative  = cumsum(co),
         x           = step_info$shift
  )

}


#' response_distributed_lag
#'
#' @description Calculate the response for a distributed lag set
#'
#' @inheritParams get_coefficients
#'
#' @return A `tibble` containing the responses from distributed lag terms
#'
#' @noRd
#'
response_distributed_lag <- function(co, step_info) {

  bl    <- step_info$basis_mat[[1]]
  co    <- get_coefficients(co, step_info)

  resp  <- bl %*% co

  tibble(coefficient = as.numeric(resp),
         cumulative  = cumsum(as.numeric(resp)),
         x           = step_info$min_knots:step_info$max_knots)

}


#' response_harmonic
#'
#' @description Calculate the response (amplitude and phase) for harmonics
#'
#' @inheritParams get_coefficients
#'
#' @return a `tibble` containing the responses from harmonic terms
#'
#' @noRd
#'
response_harmonic <- function(co, step_info) {

  co              <- get_coefficients(co, step_info)
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
#'  from `step_distributed_lag`, `step_lead_lag`, `step_harmonic` and
#'  `step_earthtide`.
#'
#' @inheritParams predict_terms
#' @param verbose Print names of steps with no response methods
#'
#' @details `step_distributed_lag` and `step_lead_lag` result
#'   in impulse response functions and `step_harmonic` and `step_earthtide`
#'   result in harmonic components (amplitude and phase for each _main_ frequency).
#'
#' @return A data.frame of impulse response functions, or harmonic components
#'   corresponding to each step.
#'
#' @export
#' @examples
#' data(transducer)
#' transducer$datetime_num <- as.numeric(transducer$datetime)
#'
#' rec_toll_rasmussen <- recipe(wl~baro + datetime_num, transducer) |>
#'    step_lead_lag(baro, lag = log_lags(100, 86400 * 2 / 120)) |>
#'    step_ns(datetime_num, deg_free = 10) |>
#'    prep()
#'
#' input_toll_rasmussen <- rec_toll_rasmussen |> bake(new_data = NULL)
#'
#' fit_toll_rasmussen <- lm(wl~., input_toll_rasmussen)
#' resp <- response(fit_toll_rasmussen,
#'                  rec_toll_rasmussen)
#' plot(value~x, resp[resp$name == 'cumulative',], type = 'l')
response <- function(fit, rec, verbose = FALSE, ...) UseMethod("response")


#' @rdname response
#' @export
response.lm <- function(fit, rec, verbose = FALSE, ...) {

  co <- coefficients(fit)

  response.numeric(co, rec, verbose = verbose,...)

}


#' @rdname response
#' @export
response.cv.glmnet <- function(fit, rec, verbose = FALSE, ...) {

  co        <- coefficients(fit)
  co_names  <- rownames(co)
  co        <- as.vector(co)
  names(co) <- co_names

  response.numeric(co, rec, verbose = verbose, ...)
}


#' @rdname response
#' @export
response.numeric <- function(fit, rec, verbose = FALSE, ...) {

  x <- NULL

  rec_steps   <- tidy2(rec)
  resp        <- vector(mode = "list", length = nrow(rec_steps))
  names(resp) <- paste0(rec_steps$type, '_', rec_steps$step_name)

  for (i in 1:nrow(rec_steps)) {

    step_info <- tidy2(rec, i)
    type      <- rec_steps$type[i]

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

      out       <- pivot_longer(out, cols = !x)
      out$type  <- type
      out$term  <- unique(step_info$terms)
      resp[[i]] <- out

    }

  }

  bind_rows(resp)
}
