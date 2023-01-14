#include "hydrorecipes.h"


//==============================================================================
//' @title
//' b_spline
//'
//' @description
//' Calculate the basis splines
//'
//' @inheritParams splines::bs
//' @param knots location of knots for the b-splines. Unlike `splines::bs` this
//' includes the boundary knots. (numeric vector)
//'
//'
//' @return the basis spline values with intercept.
//'
//' @noRd
//'
// [[Rcpp::export]]
 Eigen::MatrixXd b_spline(const Eigen::ArrayXd& x,
                          const Eigen::ArrayXd& knots,
                          size_t degree) {

   size_t n = x.size();
   size_t n_knot = knots.size();
   size_t n_cols = n_knot + degree - 1;
   size_t n_knots = (degree) * 2 + n_knot;

   size_t order = degree + 1;
   size_t k_offset = 0;
   size_t j_index = 0;

   double i1 = 0;
   double i2 = 0;

   double saved;
   double den;
   double term;

   ArrayXd knots_pad(n_knots);
   VectorXd pad(degree);

   pad.setConstant(knots(0));
   knots_pad.head(degree) = pad;
   pad.setConstant(knots(n_knot - 1));
   knots_pad.tail(degree) = pad;
   knots_pad.segment(degree, n_knot) = knots;

   MatrixXd out = MatrixXd::Zero(n, n_cols);
   VectorXi ind = which_indices(x, knots);

   // These are the indices to keep
   for (size_t i = 0; i < n; ++i) {
     out(i, ind(i)) = 1;
   }


   // the degree of the curve
   for (size_t k = 1; k <= degree; ++k) {
     k_offset = degree - k;

     // loop through each x value
     for (size_t i = 0; i < n; ++i) {
       saved = 0;

       for (size_t j = 0; j < k; ++j) {
         j_index = ind(i) + j;
         i1 = knots_pad(j_index + k_offset + 1);
         i2 = knots_pad(j_index + order);
         den  = i2 - i1;
         if(den == 0) {
           term = 0;
         } else {
           term = out(i, j_index) / den;
         }
         out(i, j_index) = saved + (i2 - x(i)) * term;
         saved = (x(i) - i1) * term;
       }
       out(i, ind[i] + k) = saved;
     }
   }

   return(out);

 }
//==============================================================================


/*** R
*/
