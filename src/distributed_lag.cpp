#include "frecipes.h"

//==============================================================================
//' @title
//' distributed_lag_thread
//'
//' @description
//' Create distributed lag terms
//'
//' @param x
//' @param bl
//' @param n_thread
//'
//' @return List of distributed lags
//'
//' @export
//'
// [[Rcpp::export]]
Eigen::MatrixXd distributed_lag_thread(const Eigen::VectorXd& x,
                                       const Eigen::MatrixXd& bl,
                                       unsigned int n_thread) {

  // result matrix
  unsigned int n_row = x.size();
  unsigned int n_col = bl.rows();
  unsigned int n_rem = bl.cols();

  Eigen::MatrixXd cb(n_row, n_col);
  cb.setConstant(NA_REAL);

  RcppThread::ThreadPool pool(n_thread);

  pool.parallelFor(0, n_row - n_rem + 1, [&] (unsigned int k) {
    cb.row(n_row - k - 1) =  bl * x.segment(k, n_rem);
  });

  pool.join();

  return cb;
}


// [[Rcpp::export]]
Eigen::VectorXd convolve_eigen(const Eigen::VectorXd& x,
                               const Eigen::RowVectorXd& y) {

  unsigned int n_row = x.size();
  unsigned int n_rem = y.size();

  Eigen::VectorXd out(n_row);
  out.setConstant(NA_REAL);

  for (unsigned int k = 0; k < n_row - n_rem + 1; ++k) {
    out(n_row - k - 1) = y.dot(x.segment(k, n_rem));
  }

  return(out);
}


//==============================================================================
//' @title
//' distributed_lag_eigen
//'
//' @description
//' Create distributed lag terms
//'
//' @param x
//' @param bl
//'
//' @return List of distributed lags
//'
//' @export
//'
// [[Rcpp::export]]
Rcpp::List distributed_lag_eigen(Eigen::Map<Eigen::VectorXd> x,
                                 Eigen::Map<Eigen::MatrixXd> bl) {

  // result matrix
  unsigned int n_row = x.size();
  unsigned int n_col = bl.rows();
  unsigned int n_rem = bl.cols();

  Rcpp::List cb(n_col);

  for (unsigned int i = 0; i < n_col; ++i){
    cb(i) = convolve_eigen(x, bl.row(i));
  };

  return cb;
}



// // [[Rcpp::export]]
// NumericVector convolve_arma(const arma::rowvec x,
//                             const arma::vec y) {
//
//   size_t n_row = x.n_elem;
//   size_t n_rem = y.n_elem;
//
//   arma::vec out(n_row);
//   out.fill(NA_REAL);
//
//   for (size_t k = 0; k < n_row - n_rem + 1; ++k) {
//     out(n_row - k - 1) = as_scalar(x.subvec(k, k+n_rem-1) * y);
//   }
//
//   return(Rcpp::NumericVector(out.begin(), out.end()));
// }
//
//
// // [[Rcpp::export]]
// Rcpp::List distributed_lag_arma(const arma::rowvec& x,
//                                 const arma::mat& bl) {
//
//   // result matrix
//   size_t n_row = x.n_elem;
//   size_t n_col = bl.n_cols;
//   size_t n_rem = bl.n_rows;
//
//   Rcpp::List cb(n_col);
//
//   for (size_t i = 0; i < n_col; ++i){
//     cb(i) = convolve_arma(x, bl.col(i));
//   };
//
//   return cb;
// }



//==============================================================================
//' @title
//' distributed_lag_thread_eigen
//'
//' @description
//' Create distributed lag terms
//'
//' @param x
//' @param bl
//' @param lag_max
//' @param n_subset
//' @param n_shift
//' @param n_thread
//'
//' @return List of distributed lags
//'
//' @export
//'
// [[Rcpp::export]]
List distributed_lag_thread_eigen(Eigen::Map<Eigen::VectorXd> x,
                                  Eigen::Map<Eigen::MatrixXd> bl,
                                  int lag_max,
                                  int n_subset,
                                  int n_shift,
                                  int n_thread) {

  // result matrix
  int n_row = x.size();
  int n_col = bl.rows();

  // added
  int n_out;
  int start;
  int end;
  int offset = 0;

  if (n_subset < 1) {
    throw std::range_error("n_subset should be 1 or greater.");
  }
  if (n_shift >= (n_subset)) {
    throw std::range_error("The absolute value of n_shift should be less than n_subset.");
  }
  if (n_shift < 0) {
    throw std::range_error("n_shift should be positive.");
  }
  if ((n_shift + lag_max) > n_row) {
    throw std::range_error("n_shift + lag_max must be less than the length of x");
  }

  n_out = get_length(n_row, n_subset);
  start = get_start(n_out, lag_max, n_subset);
  end   = get_end(n_row, n_out, lag_max, n_subset);

  Eigen::MatrixXd cb(n_col, n_out);
  cb = cb.setConstant(NA_REAL);


  if (n_subset != 1) {
    if (n_row % n_subset == 0) {
      offset = n_subset - 1;
    } else {
      offset = abs(n_row - ((n_row / n_subset) * n_subset + 1));
    }
  }

  int wh = (n_out - start - 1) * n_subset + offset - n_shift + lag_max;

  if (wh < (n_row - n_subset)) {
    start = start - 1;
  }
  offset = offset - n_shift;

  RcppThread::ThreadPool pool(n_thread);

  pool.parallelFor(n_out - end, n_out - start, [&] (size_t i) {

    int wh = (i * n_subset) + offset;

    cb.col(i) = bl * x.segment(wh, lag_max);

  });
  pool.join();

  List out;

  for (size_t i = 0; i < cb.rows(); ++i) {
    out.push_back(cb.row(i));
  }

  return(out);

}


//==============================================================================
//' @title
//' distributed_lag_list
//'
//' @description
//' Create distributed lag terms
//'
//' @param x
//' @param n_lag
//' @param max_lag
//' @param df
//' @param degree
//' @param internal_knots
//' @param boundary_knots
//' @param complete_basis
//' @param periodic
//' @param derivs
//' @param integral
//'
//' @return List of distributed lags
//'
//' @export
//'
// [[Rcpp::export]]
Rcpp::List distributed_lag_list(Eigen::Map<Eigen::VectorXd> x,
                                arma::uword n_lag,
                                arma::uword max_lag,
                                const unsigned int df,
                                const unsigned int degree,
                                const arma::vec& internal_knots,
                                const arma::vec& boundary_knots,
                                const bool complete_basis,
                                const bool periodic,
                                const unsigned int derivs,
                                const bool integral
) {

  arma::vec rng = arma::linspace(0, max_lag, max_lag + 1);
  arma::vec knots = log_lags_arma(n_lag, max_lag);
  arma::uvec one_n = { 0, n_lag - 1 };

  Rcpp::List s = b_spline_list(rng,
                               df,
                               degree,
                               knots.subvec(1, n_lag - 2),
                               knots(one_n),
                               complete_basis,
                               periodic,
                               derivs,
                               integral);


  return(convolve_list(x, s, true, true));
}

// [[Rcpp::export]]
std::list<Eigen::VectorXd> distributed_lag_list2(Eigen::Map<Eigen::VectorXd> x,
                                arma::uword n_lag,
                                arma::uword max_lag,
                                const unsigned int df,
                                const unsigned int degree,
                                const arma::vec& internal_knots,
                                const arma::vec& boundary_knots,
                                const bool complete_basis,
                                const bool periodic,
                                const unsigned int derivs,
                                const bool integral
) {

  arma::vec rng = arma::linspace(0, max_lag, max_lag + 1);
  arma::vec knots = log_lags_arma(n_lag, max_lag);
  arma::uvec one_n = { 0, n_lag - 1 };

  std::list<Eigen::VectorXd> s = b_spline_list3(rng,
                               df,
                               degree,
                               knots.subvec(1, n_lag - 2),
                               knots(one_n),
                               complete_basis,
                               periodic,
                               derivs,
                               integral);


  return(convolve_list2(x, s, true, true));
}


// [[Rcpp::export]]
Rcpp::List distributed_lag_list3(Eigen::VectorXd x,
                                 arma::uword n_lag,
                                 arma::uword max_lag,
                                 const unsigned int df,
                                 const unsigned int degree,
                                 const arma::vec& internal_knots,
                                 const arma::vec& boundary_knots,
                                 const bool complete_basis,
                                 const bool periodic,
                                 const unsigned int derivs,
                                 const bool integral
) {

  arma::vec rng = arma::linspace(0, max_lag, max_lag + 1);
  arma::vec knots = log_lags_arma(n_lag, max_lag);
  arma::uvec one_n = { 0, n_lag - 1 };

  Rcpp::List s = b_spline_list(rng,
                               df,
                               degree,
                               knots.subvec(1, n_lag - 2),
                               knots(one_n),
                               complete_basis,
                               periodic,
                               derivs,
                               integral);

  int n_x = x.size();

  if (n_x < max_lag * 10) {
    return(convolve_list(x, s, true, true));
  }

  return(convolve_overlap_save_list(x, s, 0));
}

// [[Rcpp::export]]
Rcpp::List distributed_lag_list4(Eigen::VectorXd x,
                                 Rcpp::List s,
                                 unsigned int max_lag
) {

  unsigned int n_x = x.size();

  if (n_x < max_lag * 30) {
    return(convolve_list(x, s, true, true));
  }

  return(convolve_overlap_save_list(x, s, 0));
}



/*** R

n <- 1e6
nr <- 6
m3 <- rnorm(n)
m4 <- matrix(rep(m3, nr), ncol = n)
y  <- rnorm(1e7)
bench::mark(
  dl1 <- distributed_lag_thread_eigen(y, m4, n, 1, 0, 8),
  dl2 <- distributed_lag_thread(y, m4, 8),
  dl3 <- distributed_lag_eigen(y, m4),
  # dl4 <- distributed_lag_arma(y, m4),
  check = FALSE
)


x <- rnorm(1e7)
y <- as.numeric(0:86400)
m <- matrix(rep(y, 20), ncol = 20)
l <- list(y,y,y,y,y,y,
          y,y,y,y,y,y,
          y,y,y,y,y,y,
          y,y)

n_lags <- 20
max_lag <- 86400
ll <- as.numeric(frecipes:::log_lags_arma(n_lags, max_lag))

sp <- frecipes:::b_spline_list3(y, df = 0L, degree = 3L,internal_knots = ll[2:19], boundary_knots = c(ll[1], ll[length(ll)]))
bench::mark(
  a <- frecipes:::convolve_list(x, sp, TRUE, TRUE),
  b <- frecipes:::convolve_list2(x, sp, TRUE, TRUE),
  check = FALSE
)


bench::mark(
  # frecipes:::b_spline_list(y, df = 0L, degree = 3L,internal_knots = ll[2:19], boundary_knots = c(ll[1], ll[length(ll)])),
  # frecipes:::b_spline_list2(y, df = 0L, degree = 3L,internal_knots = ll[2:19], boundary_knots = c(ll[1], ll[length(ll)])),
# ((frecipes:::distributed_lag_list(x, 20, 1e5, 0, 3, ll[2:19], c(ll[1], ll[length(ll)]), TRUE, FALSE, 0, FALSE))[[1]]),
((frecipes:::distributed_lag_list3(x,
                                   n_lags,
                                   max_lag,
                                   0,
                                   3,
                                   ll[2:(n_lags-1)],
                                   c(ll[1], ll[length(ll)]),
                                   TRUE, FALSE, 0, FALSE))[[1]]),
#a <- frecipes:::convolve_list(x, l, TRUE, TRUE),
# a <- frecipes:::convolve_list(x, l, FALSE, FALSE),
# a <- frecipes:::convolve_matrix(x, m, TRUE, TRUE),
# test(x,y,ll),
check = FALSE,
min_iterations = 1
)

a <- rnorm(nextn(1e6))
b <- rnorm(nextn(1e7))
bench::mark(fftw::FFT(a), fftw::FFT(b), check = FALSE)


# bench::mark(frecipes:::convolve_overlap_save_list(x, k, TRUE),
#             frecipes:::convolve_overlap_save_list(x, k, FALSE),
#             check = FALSE)

*/
