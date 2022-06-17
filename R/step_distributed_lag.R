#' Create a distributed lagged predictor
#'
#' `step_distributed_lag` creates a *specification* of a recipe step that
#'   will add new columns of lagged data. Lagged data will
#'   by default include NA values where the lag was induced.
#'   These can be removed with [step_naomit()], or you may
#'   specify an alternative filler value with the `default`
#'   argument.
#'
#' @inheritParams recipes::step_pca
#' @inheritParams recipes::step_center
#' @inheritParams splines::ns
#' @param basis_mat The matrix of basis kernels to convolve. This is
#'  `NULL` until computed by [prep.recipe()].
#' @param prefix A prefix for generated column names, default to "lag_".
#' @param columns A character string of variable names that will
#'  be populated (eventually) by the `terms` argument.
#' @param default Passed to `dplyr::lag`, determines what fills empty rows
#'   left by lagging (defaults to NA).
#' @details The step assumes that the data are already _in the proper sequential
#'  order_ for lagging.
#' @family row operation steps
#' @export
#' @rdname step_distributed_lag
#'
#' @examples
#' data(transducer)
#'
#' recipe(~ ., data = transducer) |>
#'   step_distributed_lag(baro, knots = log_lags(10, 86400 * 2 / 120)) |>
#'   prep() |>
#'   bake(transducer)
#'
step_distributed_lag <-
  function(recipe,
           ...,
           role = "predictor",
           trained = FALSE,
           knots = NULL,
           basis_mat = NULL,
           default = NA,
           prefix = "distributed_lag_",
           keep_original_cols = FALSE,
           columns = NULL,
           skip = FALSE,
           id = rand_id("distributed_lag")) {

    if(length(unique(knots)) < length(knots)) {
      rlang::warn("step_distributed_lag should have uniquely valued 'knots'.  Taking unique values")
      knots <- unique(knots)
    }
    if(any(knots < 0)) {
      rlang::abort("step_distributed_lag requires 'knots' argument to be greater than or equal to 0")
    }
    if(length(knots) < 2) {
      rlang::abort("step_distributed_lag requires at least two 'knots'")
    }



    add_step(
      recipe,
      step_distributed_lag_new(
        terms = enquos(...),
        role = role,
        trained = trained,
        knots = knots,
        basis_mat = basis_mat,
        default = default,
        prefix = prefix,
        keep_original_cols = keep_original_cols,
        columns = columns,
        skip = skip,
        id = id
      )
    )
  }

step_distributed_lag_new <-
  function(terms, role, trained, knots, basis_mat, default, prefix, keep_original_cols,
           columns, skip, id) {
    step(
      subclass = "distributed_lag",
      terms = terms,
      role = role,
      trained = trained,
      knots = knots,
      basis_mat = basis_mat,
      default = default,
      prefix = prefix,
      keep_original_cols = keep_original_cols,
      columns = columns,
      skip = skip,
      id = id
    )
  }


#' basis_lag
#'
#' Depending on the vector to lag and the maximum knot, distributed_lag will
#' either use an FFT (no NA and large maximum knot), or a parallel method
#' (NA, or small maximum knot).
#'
#' @inheritParams splines::ns
#'
#' @return matrix with distributed lag terms
#'
#' @importFrom splines ns
#'
#' @export
#'
basis_lag <- function(knots) {

  # generate basis functions
  max_knot <- max(knots)
  n_knots  <- length(knots)

  # generate basis lag
  splines::ns(min(knots):max_knot,
             knots = knots[-c(1, n_knots)],
             Boundary.knots = range(knots),
             intercept = TRUE)
}


#' distributed_lag
#'
#' Depending on the vector to lag and the maximum knot, distributed_lag will
#' either use an FFT (no NA and large maximum knot), or a parallel method
#' (NA, or small maximum knot).
#'
#' @inheritParams splines::ns
#' @param basis_mat lag matrix for convolution
#'
#' @return matrix with distributed lag terms
#'
#' @importFrom splines ns
#'
#' @export
distributed_lag <- function(x, basis_mat, knots) {
  # x_len    <- length(as.numeric(x))
  max_knot <- max(knots)

  # convolution - fft for large number of lags, otherwise use parallel version
  if(any(is.na(x)) | max_knot < 5000) {
    dist_lag_mat <- distributed_lag_parallel(rev(x),
                                             t(as.matrix(basis_mat)),
                                             max(knots) - min(knots),
                                             n_subset = 1,
                                             n_shift = 0
    )
  } else {
    dist_lag_mat <- cross_basis_fft(as.matrix(x), basis_mat)
  }

  return(dist_lag_mat)
}



#' convolve_fft
#'
#' @param x numeric vector
#' @param y convolution kernel
#'
#' @return numeric vector result of convolution
#' @export
#'
convolve_fft <- function(x, y)
{
  n_in <- length(x)
  ny   <- length(y)
  n1   <- ny - 1
  x    <- c(rep.int(0, n1), x)
  n    <- length(y <- c(y, rep.int(0, n_in - 1)))
  x    <- fftw::IFFT(fftw::FFT(x) * Conj(fftw::FFT(y)), scale = FALSE) / n
  x[1:n1] <- NA_real_
  return(Re(x)[1:n_in])
}



#' cross_basis_fft
#'
#' @param basis_var numeric vector
#' @param basis_mat lagging matrix
#'
#' @return numeric vector result of convolution
#' @export
#'
cross_basis_fft <- function(basis_var, basis_mat)
{

  mat <- matrix(NA_real_,
                nrow = nrow(basis_var),
                ncol = ncol(basis_var) * ncol(basis_mat))


  for(v in seq(length = ncol(basis_var))) {
    for(l in seq(length = ncol(basis_mat))) {

      mat[, ncol(basis_mat) * (v-1)+l] <-
        convolve_fft(basis_var[,v], rev(basis_mat[,l]))

    }
  }

  return(mat)

}




#' @export
prep.step_distributed_lag <- function(x, training, info = NULL, ...) {

  col_names <- recipes_eval_select(x$terms, training, info)

  x_len <- nrow(training)

  if(max(x$knots) > x_len) {
    rlang::abort('The maximum knot cannot be larger than the number of elements in x')
  }

  basis_mat <- basis_lag(x$knots)

  step_distributed_lag_new(
    terms = x$terms,
    role = x$role,
    trained = TRUE,
    knots = x$knots,
    basis_mat = basis_mat,
    default = x$default,
    prefix = x$prefix,
    keep_original_cols = x$keep_original_cols,
    columns = recipes_eval_select(x$terms, training, info),
    skip = x$skip,
    id = x$id
  )
}


#' @export
bake.step_distributed_lag <- function(object, new_data, ...) {

  for(i in seq_along(object$columns)) {

    dl <- distributed_lag(new_data[[object$columns[i]]],
                          object$basis_mat,
                          object$knots)


    colnames(dl) <- paste0(object$prefix, object$columns[i], '_',
                           object$knots)


    new_data <- bind_cols(new_data, as_tibble(dl))

  }

  keep_original_cols <- get_keep_original_cols(object)
  if (!keep_original_cols) {
    new_data <-
      new_data[, !(colnames(new_data) %in% object$columns), drop = FALSE]
  }

  new_data

}

print.step_distributed_lag <-
  function(x, width = max(20, options()$width - 30), ...) {
    cat("Distributed lag model with ",  sep = "")
    printer(x$columns, x$terms, x$trained, width = width)
    invisible(x)
  }


#' @export
tidy.step_distributed_lag <- function(x, ...) {

  tidy2.step_distributed_lag(x, ...)
}

#' @export
tidy2.step_distributed_lag <- function(x, ...) {

  if (is_trained(x)) {
    res <-
      tibble(terms = rep(x$columns, each = length(x$knots)),
             knots = rep(x$knots, times = length(x$columns)))
    res$key <- paste0(x$prefix, res$terms, '_', res$knots)

  } else {
    term_names <- sel2char(x$terms)

    res <- tibble(terms = rep(term_names, each = length(x$knots),
                              knots = rep(x$knots, times = length(term_names))))
    res$key <- paste0(x$prefix, res$terms, '_',
                      rep(x$knots, times = length(term_names)))

  }

  res$id <- x$id
  res$step_name <- 'step_distributed_lag'

  res
}
