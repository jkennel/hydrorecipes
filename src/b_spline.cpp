#include "hydrorecipes.h"


// #define EIGEN_FFTW_DEFAULT
//
// #include <RcppArmadillo.h>
// #include <RcppEigen.h>
// #include <RcppParallel.h>
// #include <RcppThread.h>
// // #include <fftw3.h>
// // #include <unsupported/Eigen/FFT>
//
//
// using namespace RcppParallel;
// using namespace arma;
// using namespace Rcpp;
//
// // [[Rcpp::depends(RcppArmadillo)]]
// // [[Rcpp::depends(RcppEigen)]]
// // [[Rcpp::depends(RcppThread)]]
//
//
// using Eigen::MatrixXd;                  // variable size matrix, double precision
// using Eigen::VectorXd;
// using Eigen::VectorXi;
// using Eigen::ArrayXd;
//
//
// // // [[Rcpp::export]]
// // IntegerVector fi(const NumericVector& x,
// //                  const NumericVector& vec,
// //                  bool rightmost_closed,
// //                  bool all_inside,
// //                  bool left_open) {
// //
// //   Rcpp::Function f("findInterval");
// //
// //   const IntegerVector out = f(x, vec, rightmost_closed, all_inside, left_open);
// //
// //   return(out);
// //
// // }
// //
// // // [[Rcpp::export]]
// // Eigen::VectorXi fi2(const Eigen::VectorXd& x,
// //                     const Eigen::VectorXd& vec,
// //                     bool rightmost_closed,
// //                     bool all_inside,
// //                     bool left_open) {
// //
// //   Rcpp::Function f("findInterval");
// //
// //   IntegerVector tmp = f(wrap(x), wrap(vec),
// //                         rightmost_closed, all_inside,
// //                         left_open);
// //   Eigen::VectorXi out = as<Eigen::VectorXi>(tmp);
// //
// //   return(out);
// //
// // }
// // [[Rcpp::export]]
// Eigen::VectorXi which_indices(const Eigen::VectorXd& x,
//                               const Eigen::VectorXd& knots) {
//
//   size_t n_x = x.size();
//   size_t n_below = 0;
//   size_t n_old = 0;
//   size_t n_i = 0;
//   size_t n = 0;
//   size_t n_knots = knots.size();
//   Eigen::VectorXi ind = Eigen::VectorXi::Constant(n_x, n_knots - 2);
//
//
//   for (size_t i = 0; i < n_knots; ++i) {
//
//     n_i = (x.array() <= knots(i)).count();
//     n_below = n_i - n_old;
//
//     if (n_below > 0) {
//       if (i == 0) {
//         ind.segment(n, n_below).setConstant(i);
//       } else {
//         ind.segment(n, n_below).setConstant(i-1);
//       }
//     }
//
//     n_old = n_i;
//     n += n_below;
//
//   }
//
//   return(ind);
// }
// //==============================================================================
// //' @title
// //' b_spline
// //'
// //' @description
// //' Calculate the basis splines
// //'
// //' @inheritParams splines::bs
// //' @param knots location of knots for the b-splines. Unlike `splines::bs` this
// //' includes the boundary knots. (numeric vector)
// //'
// //'
// //' @return the basis spline values with intercept.
// //'
// //' @noRd
// //'
// // [[Rcpp::export]]
// std::list<Eigen::VectorXd> b_spline(Eigen::Map<Eigen::VectorXd> x,
//                           Eigen::Map<Eigen::VectorXd> knots,
//                           size_t degree) {
//
//    size_t n = x.size();
//    size_t n_knot = knots.size();
//    size_t n_cols = n_knot + degree - 1;
//    size_t n_knots = (degree) * 2 + n_knot;
//
//    size_t order = degree + 1;
//    size_t k_offset = 0;
//    size_t j_index = 0;
//
//    double i1 = 0;
//    double i2 = 0;
//
//    double saved;
//    double den;
//    double term;
//
//    ArrayXd knots_pad(n_knots);
//    VectorXd pad(degree);
//
//    pad.setConstant(knots(0));
//    knots_pad.head(degree) = pad;
//    pad.setConstant(knots(n_knot - 1));
//    knots_pad.tail(degree) = pad;
//    knots_pad.segment(degree, n_knot) = knots;
//
//    MatrixXd out = MatrixXd::Zero(n, n_cols);
//    VectorXi ind = which_indices(x, knots);
//
//    // These are the indices to keep
//    for (size_t i = 0; i < n; ++i) {
//      out(i, ind(i)) = 1;
//    }
//
//    // the degree of the curve
//    for (size_t k = 1; k <= degree; ++k) {
//      k_offset = degree - k;
//
//      // loop through each x value
//      for (size_t i = 0; i < n; ++i) {
//        saved = 0;
//
//        for (size_t j = 0; j < k; ++j) {
//          j_index = ind(i) + j;
//          i1 = knots_pad(j_index + k_offset + 1);
//          i2 = knots_pad(j_index + order);
//          if (i1 == i2) {
//            term = 0;
//          } else {
//            den  = i2 - i1;
//            term = out(i, j_index) / den;
//            out(i, j_index) = saved + (i2 - x(i)) * term;
//            saved = (x(i) - i1) * term;
//          }
//        }
//        out(i, ind(i) + k) = saved;
//      }
//    }
//
//    std::list<Eigen::VectorXd> o;
//
//    for (size_t i = 0; i < out.cols(); ++i) {
//      o.push_back(out.col(i));
//    }
//
//    return(o);
//
//  }
//
//
// // [[Rcpp::export]]
// NumericMatrix b_spline2(const NumericVector& x,
//                         const NumericVector& knots,
//                         size_t degree) {
//
//   size_t n = x.size();
//   size_t n_knot = knots.size();
//   size_t n_cols = n_knot + degree - 1;
//   size_t n_knots = (degree) * 2 + n_knot;
//
//   size_t order = degree + 1;
//   size_t k_offset = 0;
//   size_t j_index = 0;
//
//   double i1 = 0;
//   double i2 = 0;
//
//   double saved;
//   double den;
//   double term;
//
//   NumericVector knots_pad(n_knots);
//   NumericVector pad(degree);
//
//   pad.fill(knots(0));
//   knots_pad[Range(0, degree - 1)] = pad;
//   pad.fill(knots(n_knot - 1));
//   knots_pad[Range(n_knots - degree, n_knots-1)] = pad;
//   knots_pad[Range(degree, degree + n_knot - 1)] = knots;
//
//   NumericMatrix out(n,n_cols);
//   IntegerVector ind = fi(x, knots, true, true, true)-1;
//
//   // These are the indices to keep
//   for (size_t i = 0; i < n; ++i) {
//     out(i, ind(i)) = 1;
//   }
//
//
//   // the degree of the curve
//   for (size_t k = 1; k <= degree; ++k) {
//     k_offset = degree - k;
//
//     // loop through each x value
//     for (size_t i = 0; i < n; ++i) {
//       saved = 0;
//
//       for (size_t j = 0; j < k; ++j) {
//         j_index = ind(i) + j;
//         i1 = knots_pad(j_index + k_offset + 1);
//         i2 = knots_pad(j_index + order);
//         den  = i2 - i1;
//         if (den == 0) {
//           term = 0;
//         } else {
//           term = out(i, j_index) / den;
//         }
//         out(i, j_index) = saved + (i2 - x(i)) * term;
//         saved = (x(i) - i1) * term;
//       }
//       out(i, ind(i) + k) = saved;
//     }
//   }
//
//   return(out);
//
// }




/*** R
# n <- 1e7
# x <- as.numeric(1:n)
# knots <- c(10, 100, 1000, 100000)
#
# bench::mark(
#   # splines::bs(x, knots = knots),
#   a <- splines2::bSpline(x, knots = knots, intercept = TRUE),
#   b <- b_spline(x, c(1, knots, n), 3),
#   # b_spline2(x, c(1, knots, n), 3),
#   check = FALSE,
#   iterations = 5
# )


# bench::mark(
#   (fi2(x, c(knots), TRUE,TRUE,TRUE))-1L,
#   (fi(x, c(knots), TRUE,TRUE,TRUE))-1L,
#   (findInterval(x, c(knots), TRUE,TRUE,TRUE))-1L,
#   (which_indices(x, c(knots))),
#   check = FALSE
# )
# all.equal((which_indices(x, c(knots))),
#           (fi2(x, c(knots), TRUE,TRUE,TRUE)-1)
# )

*/
