#include "frecipes.h"

//==============================================================================
//' @title
//' harmonic_list
//'
//' @description
//' Create sin and cosine terms for harmonic analysis
//'
//' @param time numeric vector of times
//' @param frequency numeric vector of frequencies
//' @param start time the cycle starts
//' @param cycle_size size of the cycle in number of measurements
//'
//' @return List of cosines and sines
//'
//' @export
//'
// [[Rcpp::export]]
Rcpp::List harmonic_list(const Rcpp::NumericVector& time,
                         const Rcpp::NumericVector& frequency,
                         const double start,
                         const double cycle_size) {


  const Rcpp::NumericVector m = (M_2PI / cycle_size) * (time - start);

  Rcpp::List out;

  for (auto &f : frequency) {
    out.push_back(cos(m * f));
    out.push_back(sin(m * f));
  }

  return(out);
}




// // [[Rcpp::export]]
// List harmonic_list_2(Eigen::VectorXd& time,
//                    Eigen::VectorXd& frequency,
//                    const double start,
//                    const double cycle_size) {
//
//
//   Eigen::VectorXd m = (M_2PI / cycle_size) * (time.array() - start);
//   unsigned int n = frequency.size();
//
//   List out(n * 2);
//   unsigned int j = 0;
//
//   for (unsigned int i = 0; i < n; ++i) {
//     out[j] = sin(m.array() * frequency(i));
//     j += 1;
//     out[j] = cos(m.array() * frequency(i));
//     j += 1;
//   }
//
//   return(out);
// }



// // [[Rcpp::export]]
// arma::fmat harmonic_float(const arma::fcolvec& time,
//                           const arma::frowvec& frequency,
//                           float start,
//                           float cycle_size) {
//
//   // if (start == NA_REAL) {
//   //   start = time[0];
//   // }
//
//   arma::fmat m1 = ((M_2PI / cycle_size) * ((time - time[0]))) * frequency;
//
//   return arma::join_rows(arma::sin(m1), arma::cos(m1));
//
// }
//
// //==============================================================================
// //' harmonic
//  //'
//  //' calculate sin and cos curves from POSIXct times (serial)
//  //'
//  //' @param time \code{numeric vector} times to calculate sin and cos at
//  //' @param frequency \code{numeric vector} frequencies for sin and cos
//  //' @param start \code{double} the starting point
//  //' @param cycle_size \code{double} the size of one cycle (period).
//  //'
//  //' @return sin and cos curves
//  //'
//  //' @export
//  //'
//  // [[Rcpp::export]]
//  arma::mat harmonic(const arma::colvec time,
//                     const arma::rowvec frequency,
//                     double start,
//                     double cycle_size) {
//
//    //   if (start == NA_REAL) {
//    //     start = time[0];
//    //   }
//
//    arma::mat m = ((M_2PI / cycle_size) * ((time - time[0]))) * frequency;
//
//    return arma::join_rows(arma::sin(m), arma::cos(m));
//
//  }
//
//
// // [[Rcpp::export]]
// arma::field<arma::colvec> harmonic_field(const arma::colvec time,
//                                          const arma::colvec frequency,
//                                          double start,
//                                          double cycle_size) {
//
//   // if (start == NA_REAL) {
//   //   start = time[0];
//   // }
//
//   arma::vec m = (M_2PI / cycle_size) * (time - time[0]);
//   arma::field<arma::colvec> out(2 * frequency.n_elem);
//
//   for (size_t i = 0; i < frequency.n_elem; ++i) {
//     out[i] = arma::sin(m * frequency[i]);
//     out[i + frequency.n_elem] = arma::cos(m * frequency[i]);
//   }
//
//   return out;
//
// }
//
//
//
//
// // [[Rcpp::export]]
// List harmonic_list2(const arma::vec& time,
//                     const arma::vec& frequency,
//                     double start,
//                     double cycle_size) {
//
//
//   const arma::vec m = (M_2PI / cycle_size) * (time-start);
//   List out;
//
//   for (size_t i = 0; i < frequency.size(); ++i) {
//     out.push_back(sin(m * frequency(i)));
//     out.push_back(cos(m * frequency(i)));
//   }
//
//   return(out);
// }
//
// // [[Rcpp::export]]
// arma::vec harmonic_interp(arma::vec time) {
//   arma::vec phase = arma::regspace<vec>(-2*M_PI, M_PI/1000, 2*M_PI);
//   arma::vec res = cos(phase);
//   arma::vec out(time.n_elem);
//
//   arma::interp1(phase, res, time, out, "*linear");
//   return out;
// }
//
// // // [[Rcpp::export]]
// // Eigen::ArrayXd harmonic_eigen(Eigen::ArrayXd& time) {
// //   return time.cos();
// // }
// //
// // [[Rcpp::export]]
// NumericVector harmonic_nv(NumericVector& time) {
//   return cos(time);
// }
//
// // [[Rcpp::export]]
// Eigen::MatrixXd harmonic_double(Eigen::Map<Eigen::VectorXd>  x,
//                                 const Eigen::RowVectorXd& frequency,
//                                 double cycle_size) {
//
//   // the different cycles
//   Eigen::VectorXd cycles = (M_2PI / cycle_size) * (x.array() - x(0));
//
//   size_t n_freq = frequency.size();
//   Eigen::MatrixXd out(x.size(), n_freq * 2);
//   Eigen::ArrayXXd m1 = (cycles * frequency);
//
//   // sin curves are first columns and cos curves are the next set of columns
//   out << m1.sin(), m1.cos();
//
//   return out;
//
// }


/*** R

n <- 1e6L
time <- sort(rnorm(n))
frequency <- c(1,2,3,4,5,6,7)
cycle_size <- 86400
start <- 0

sincos <- function(time, frequency, start, cycle_size) {

  l <- list(time)
  for (i in seq_along(frequency)) {
    m <- (2.0 * pi / cycle_size) * (time - start) * frequency[i]
    add_vars(l, list(sin(m), cos(m)))
  }

}


bench::mark(
  h0 <- frecipes:::harmonic_list(time, frequency, 0, 86400),
  # h1 <- sincos(time, frequency, 0, 86400),
  # h1 <- frecipes:::harmonic_list_2(t, vec, 0, 86400),
  # h1 <- frecipes:::harmonic_std_list(t, vec, 0, 86400),
  check = FALSE
)

*/
