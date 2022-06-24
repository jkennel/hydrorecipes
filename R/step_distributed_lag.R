#' Create a distributed lagged predictor
#'
#' `step_distributed_lag` creates a *specification* of a recipe step that
#'   will add new basis lag columns. The new data will
#'   include NA values up to the maximum lag. These can be removed
#'   with [recipes::step_naomit()]. The inspiration for this step comes from the
#'   [dlnm package](https://CRAN.R-project.org/package=dlnm). For large datasets
#'   with large maximum time lags, convolution is
#'   done in the frequency domain for efficiency. Samples should be ordered and
#'   have regular spacing (i.e. regular time series, regular spatial sampling).
#'
#' @inheritParams recipes::step_pca
#' @inheritParams recipes::step_center
#' @param knots An integer vector of breakpoints to define the spline. These should
#'  include the `Boundary.knots`. See [splines](https://CRAN.R-project.org/package=splines)
#'  for more info.
#' @param basis_mat The matrix of basis kernels to convolve. This is
#'  `NULL` until computed by [prep.recipe()]. This can also be specified as an
#'  object generated from the [splines](https://CRAN.R-project.org/package=splines) or [splines2](https://CRAN.R-project.org/package=splines2) packages having attributes
#'  for `knots` and `Boundary.knots`. If specified like this `knots` will be
#'  obtained from the `basis_mat` and not from the `knots` parameter.
#' @param spline_fun Function used for calculating `basis_mat`.  This should
#'  return an object having `knots` and `Boundary.knots` attributes.
#' @param options The arguments to pass to `spline_fun`.
#' @param prefix A prefix for generated column names, default to "distributed_lag_".
#' @family row operation steps
#' @export
#' @rdname step_distributed_lag
#'
#' @details This step assumes that the data are already _in the proper sequential
#'  order_ for lagging. The input should be sampled at a regular interval
#'  (time, space, etc.). When the recipe is baked a set of vectors resulting from
#'  the convolution of a vector and a basis matrix is returned. Distributed lags
#'  can be used to model a delayed response to a input
#'  in a flexible manner with fewer regressor terms. The method achieves this by
#'  convolving a input stress with a basis lag matrix (commonly spline function)
#'  which leads to a set of regressors with fewer terms but still capable of
#'  describing both fast and slow responses.
#'
#'
#' @return An updated version of recipe with the new step added to the sequence
#' of any existing operations.
#'
#' @references
#' Almon, S (1965). The Distributed Lag Between Capital Appropriations and
#'  Expenditures. Econometrica 33(1), 178.
#'
#' Gasparrini A. Distributed lag linear and non-linear models in R: the package
#'  dlnm. Journal of Statistical Software. 2011; 43(8):1-20.
#'  https://doi.org/10.18637/jss.v043.i08
#'
#' @examples
#' data(wipp30)
#'
#' rec_base <- recipe(wl~baro, data = wipp30)
#'
#' # default uses splines::ns
#' rec <- rec_base |>
#'   step_distributed_lag(baro,
#'                        knots = log_lags(4, 72)) |>
#'   prep()
#'
#' # use different spline function
#' rec <- rec_base |>
#'   step_distributed_lag(baro,
#'                        spline_fun = splines::bs,
#'                        options = list(intercept = TRUE,
#'                                       degree = 4L),
#'                        knots = log_lags(4, 72)) |>
#'   prep()
#'
#' # specify basis_mat
#' basis_mat <- splines2::mSpline(0:72, knots = c(3,16))
#' rec <- rec_base |>
#'   step_distributed_lag(baro,
#'                        basis_mat = basis_mat) |>
#'   prep()
#'
#' @seealso [step_lead_lag()] [recipes::step_lag()]
step_distributed_lag <-
  function(recipe,
           ...,
           role = "predictor",
           trained = FALSE,
           knots = NULL,
           basis_mat = NULL,
           spline_fun = splines::ns,
           options = list(intercept = TRUE),
           prefix = "distributed_lag_",
           keep_original_cols = FALSE,
           columns = NULL,
           skip = FALSE,
           id = rand_id("distributed_lag")) {

    if(length(unique(knots)) < length(knots)) {
      # rlang::warn("step_distributed_lag should have uniquely valued 'knots'.  Taking unique values")
      knots <- unique(knots)
    }
    if(any(knots < 0)) {
      rlang::abort("step_distributed_lag requires 'knots' argument to be greater than or equal to 0")
    }
    if(length(knots) < 2) {
      if(is.null(basis_mat)) {
        rlang::abort("step_distributed_lag requires at least two 'knots'")
      }
    }
    if(!all(as.integer(knots) == knots)) {
      rlang::warn("step_distributed_lag should have integer valued 'knots'")
    }



    add_step(
      recipe,
      step_distributed_lag_new(
        terms = enquos(...),
        role = role,
        trained = trained,
        knots = knots,
        basis_mat = basis_mat,
        spline_fun = spline_fun,
        options = options,
        prefix = prefix,
        keep_original_cols = keep_original_cols,
        columns = columns,
        skip = skip,
        id = id
      )
    )
  }

step_distributed_lag_new <-
  function(terms, role, trained, knots, basis_mat, spline_fun, options, prefix,
           keep_original_cols, columns, skip, id) {
    step(
      subclass = "distributed_lag",
      terms = terms,
      role = role,
      trained = trained,
      knots = knots,
      basis_mat = basis_mat,
      spline_fun = spline_fun,
      options = options,
      prefix = prefix,
      keep_original_cols = keep_original_cols,
      columns = columns,
      skip = skip,
      id = id
    )
  }


#' basis_lag
#'
#' Create a spline basis matrix. This function passes a list of knots and
#' options to a spline function.  Common spline functions used come from the
#' `splines` and `splines2` packages.
#'
#' @inheritParams splines::ns
#' @inheritParams step_distributed_lag
#' @param options list of named variables to pass to spline_fun
#'
#' @return \code{numeric matrix} with distributed lag terms
#'
#' @importFrom splines ns
#'
#' @noRd
basis_lag <- function(knots,
                      spline_fun = splines::ns,
                      options = list(intercept = TRUE)) {

  # generate basis functions
  max_knot <- max(knots)
  n_knots  <- length(knots)

  # generate basis lag
  do.call(spline_fun,
          c(list(x = min(knots):max(knots)),
                 list(knots = knots[-c(1, n_knots)]),
                 options))

}


#' distributed_lag
#'
#' Convolves a vector or matrix with a basis matrix (a cross basis matrix).
#' Depending on the input and the maximum knot, distributed_lag will
#' either use an FFT (no NA values and large maximum knot), or a parallel method
#' (NA values present, or small maximum knot).
#'
#' @inheritParams splines::ns
#' @param basis_mat \code{numeric matrix} spline basis matrix for convolution
#'
#' @return matrix with distributed lag terms
#'
#' @importFrom splines ns
#'
#' @noRd
distributed_lag <- function(x, basis_mat, knots) {

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
#' Do an FFT convolution of two vectors using the fftw package.
#'
#' @param x numeric vector to convolve with y
#' @param y numeric vector convolution kernel
#'
#' @return numeric vector that is the result of convolution
#'
#' @noRd
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
#' Convolve each column of a input matrix (`basis_var`) with a matrix of basis
#' splines (`basis_mat`)
#'
#' @param basis_var numeric matrix of input vectors
#' @param basis_mat numeric matrix of basis splines
#'
#' @return numeric vector result of convolution
#'
#' @noRd
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



  # check if basis_mat is provided
  if(is.null(x$basis_mat)) {
    basis_mat <- basis_lag(x$knots, x$spline_fun, x$options)

    # use provided basis_mat
  } else {
    if(length(intersect(c('knots', 'Boundary.knots'),
                        names(attributes(x$basis_mat)))) == 2) {
      x$knots <- sort(unique(c(attr(x$basis_mat, 'knots'),
                               attr(x$basis_mat, 'Boundary.knots'))))
      # rlang::warn('step_distributed_lag: knots are determined from basis_mat')
    } else {
      if(nrow(x$basis_mat) != (diff(range(x$knots)) + 1)) {
        rlang::abort('step_distributed_lag: basis_mat number of rows must equal
                     the range of knots plus 1')
      }
    }
    basis_mat <- x$basis_mat
  }

  if(max(x$knots) > x_len) {
    rlang::abort('The maximum knot cannot be larger than the number of elements in x')
  }

  step_distributed_lag_new(
    terms = x$terms,
    role = x$role,
    trained = TRUE,
    knots = x$knots,
    basis_mat = basis_mat,
    spline_fun = x$spline_fun,
    options = x$options,
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
                           1:ncol(object$basis_mat))


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


# We want to store the basis_matrix for
#' @rdname tidy2.recipe
#' @export
tidy2.step_distributed_lag <- function(x, ...) {

  if (is_trained(x)) {
    term_names = x$columns
  } else {
    term_names <- sel2char(x$terms)
  }

  res <- tibble(terms = term_names,
                key = paste0(x$prefix, term_names),
                min_knots = min(x$knots),
                max_knots = max(x$knots),
                spline_fun = list(x$spline_fun),
                basis_mat = list(x$basis_mat))

  res$id <- x$id
  res$step_name <- 'step_distributed_lag'

  res
}
