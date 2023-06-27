#include "hydrorecipes.h"

//==============================================================================
//' @title
//' harmonic_double
//'
//' @description
//' This function creates a matrix of sin and cosine curves at specified
//' frequencies with a defined cycle size.
//'
//' @inheritParams step_harmonic
//' @param x the sample spacings (numeric vector)
//' @param frequency frequencies of the sin and cosine curves (numeric vector)
//' @param cycle_size the cycle size in terms of x (numeric)
//'
//' @return matrix of sin and cosine curves.  The sin curves are the first set
//' of columns and the cosine curves are the second set.
//'
//' @noRd
//'
// [[Rcpp::export]]
Eigen::MatrixXd harmonic_double(const Eigen::VectorXd& x,
                                const Eigen::RowVectorXd& frequency,
                                double cycle_size) {

  // the different cycles
  VectorXd cycles = (M_2PI / cycle_size) * (x.array() - x(0));

  size_t n_freq = frequency.size();
  MatrixXd out(x.size(), n_freq * 2);
  ArrayXXd m1 = (cycles * frequency);

  // sin curves are first columns and cos curves are the next set of columns
  out << m1.sin(), m1.cos();

  return out;

}
//==============================================================================
