#include "frecipes.h"

//==============================================================================
//' @title
//' harmonic_list
//'
//' @description
//' Create sin and cosine terms for harmonic analysis
//'
//' @param time vector of times
//' @param frequency
//' @param start
//' @param cycle_size size of the cycle
//'
//' @return List of sines and cosines
//'
//' @export
//'
//'
// [[Rcpp::export]]
List harmonic_list(const NumericVector& time,
                   const NumericVector& frequency,
                   const double start,
                   const double cycle_size) {


  NumericVector m = clone(time);
  size_t n = frequency.size();

  m = (M_2PI / cycle_size) * (time - start);

  List out;

  for(size_t i = 0; i < frequency.size(); ++i) {
    out.push_back(sin(m * frequency(i)));
    out.push_back(cos(m * frequency(i)));
  }

  return(out);
}


// // [[Rcpp::export]]
// arma::fmat harmonic_float(const arma::fcolvec& time,
//                           const arma::frowvec& frequency,
//                           float start,
//                           float cycle_size) {
//
//   // if(start == NA_REAL) {
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
//    //   if(start == NA_REAL) {
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
//   // if(start == NA_REAL) {
//   //   start = time[0];
//   // }
//
//   arma::vec m = (M_2PI / cycle_size) * (time - time[0]);
//   arma::field<arma::colvec> out(2 * frequency.n_elem);
//
//   for(size_t i = 0; i < frequency.n_elem; ++i) {
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
//   for(size_t i = 0; i < frequency.size(); ++i) {
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

n <- 10000000L
t <- sort(rnorm(n))
vec <- c(1,2,3,4,5,6)

bench::mark(
  h0 <- harmonic_list(t, vec, 0, 86400),
  iterations = 3
)

*/
