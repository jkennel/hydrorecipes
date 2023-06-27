#include "hydrorecipes.h"

//==============================================================================
//' @title
//' index_from_i_j
//'
//' @description
//' Get the column number for the cross-spectra matrix.
//'
//' @param i row number (integer)
//' @param i column number (integer)
//' @param n_col the number of columns i.e. the number of covariates. (integer)
//'
//' @return integer value containing the column number of the cross-spectra.
//'
//' @noRd
//'
// [[Rcpp::export]]
size_t index_from_i_j(size_t i, size_t j, size_t n_col) {

  if(i >= n_col) {
    Rcpp::stop("index_from_i_j: i cannot be greater than or equal to n_col");
  }
  if(j >= n_col) {
    Rcpp::stop("index_from_i_j: j cannot be greater than or equal to n_col");
  }

  return(i * n_col + j);
  }
//==============================================================================

// [[Rcpp::export]]
size_t get_column_number(size_t n) {
  return((-1 + sqrt(1 + 4 * 2 * n)) / 2);
}
//==============================================================================
//' @title
//' index_from_j_i
//'
//' @description
//' Get the column number for the cross-spectra matrix.
//'
//' @inheritParams index_from_i_j
//'
//' @return integer value containing the column number of the cross-spectra
//' conjugate term.
//'
//' @noRd
//'
// [[Rcpp::export]]
size_t index_from_j_i(size_t i, size_t j, size_t n_col) {
  return(index_from_i_j(j, i, n_col));
}
//==============================================================================


//==============================================================================
//' @title
//' next_n_eigen
//'
//' @description
//' Get the size of the nearest fast length for an FFT. This just uses the base
//' function right now.
//'
//' @inheritParams nextn
//'
//' @return 'nice' integer value for the length of the padded FFT.
//'
//' @noRd
//'
// [[Rcpp::export]]
size_t next_n_eigen(size_t n) {

  // Get efficient padding
  Rcpp::Environment base("package:stats");
  Rcpp::Function nextn_r = base["nextn"];
  size_t n_new = *INTEGER(nextn_r(n));

  return(n_new);
}
//==============================================================================


//==============================================================================
//' @title
//' pad_vector
//'
//' @description
//' Pad a vector with zeros to a desired length.
//'
//' @param x initial vector to pad (numeric vector)
//' @param x_old initial length of vector (integer)
//' @param x_new desired length of padded vector (integer)
//'
//' @return 'nice' integer value for the length of the padded FFT.
//'
//' @noRd
//'
// [[Rcpp::export]]
Eigen::VectorXd pad_vector(Eigen::VectorXd x, size_t n_old, size_t n_new) {

  if (n_new < n_old) {
    Rcpp::stop("pad_vector: n_new cannot be smaller than n_old");
  } else if (n_new == n_old) {
    return(x);
  }

  x.conservativeResize(n_new);
  x.tail(n_new - n_old).setZero();

  return(x);
}
//==============================================================================


//==============================================================================
//' @title
//' detrend_matrix
//'
//' @description
//' Linearly detrend the columns of a matrix. This is translated from spec.pgram
//'
//' @param x the matrix that holds multiple series (numeric matrix)
//'
//' @return columns of a matrix that have been linearly detrended.
//'
//' @noRd
//'
// [[Rcpp::export]]
Eigen::MatrixXd detrend_matrix(const Eigen::MatrixXd& x) {

  size_t n_row = x.rows();
  size_t n_col = x.cols();

  double ends = ((double)n_row - 1.0) / 2.0;
  double scale = n_row * (n_row * n_row - 1.0) / 12.0;

  VectorXd x_abscissa = VectorXd::LinSpaced(n_row, -ends, ends);

  MatrixXd out = x;

  for (size_t i = 0; i < n_col; ++i) {

    out.col(i) = out.col(i).array() - out.col(i).array().mean() -
      (out.col(i).array() * x_abscissa.array()).sum() *
      x_abscissa.array() / scale;

  }

  return(out);
}
//==============================================================================


//==============================================================================
//' @title
//' demean_matrix
//'
//' @description
//' Remove the mean from each column of a matrix.
//'
//' @inheritParams detrend_matrix
//'
//' @return columns of a matrix with the means removed.
//'
//' @noRd
//'
// [[Rcpp::export]]
Eigen::MatrixXd demean_matrix(const Eigen::MatrixXd& x) {

  size_t n_row = x.rows();
  size_t n_col = x.cols();

  MatrixXd out(n_row, n_col);

  for (size_t i = 0; i < n_col; ++i) {
    out.col(i) = x.col(i).array() - x.col(i).mean();
  }

  return(out);
}
//==============================================================================


//==============================================================================
//' @title
//' detrend_and_demean_matrix
//'
//' @description
//' Remove the trend and mean from each column of a matrix.
//'
//' @inheritParams detrend_matrix
//' @param detrend should the trend be removed from each column (boolean)
//' @param demean should the mean be removed from each column (boolean)
//'
//' @return columns of a matrix with the means and/or trends removed.
//'
//' @noRd
//'
// [[Rcpp::export]]
Eigen::MatrixXd detrend_and_demean_matrix(const Eigen::MatrixXd& x,
                                          bool detrend,
                                          bool demean) {

  if (detrend) {
    return(detrend_matrix(x));
  } else if (demean) {
    return(demean_matrix(x));
  }

  return(x);
}
//==============================================================================




// //==============================================================================
// //' @title
// //' modified_daniell
// //'
// //' @description
// //' Create a modified daniell kernel using FFT. Adapted from `spec.pgram`.
// //'
// //' @inheritParams spec.pgram
// //'
// //' @param n length of the kernel (integer)
// //'
// //' @return modified Daniell kernel.
// //'
// //' @noRd
// //'
// // [[Rcpp::export]]
// Eigen::VectorXd modified_daniell(Eigen::VectorXi spans, size_t n) {
//
//   Eigen::VectorXd out(n);
//   out.setZero();
//
//   size_t n_kernel = spans.sum();
//   size_t n_spans = spans.size();
//   size_t m, kernel_length;
//
//   Eigen::VectorXd kernel(n_kernel*2 + 1);
//   kernel.setZero();
//   Eigen::VectorXd kernel_run(n_kernel*2 + 1);
//   kernel_run.setZero();
//
//   m = spans(0);
//   kernel_length = (m + 1);
//
//   kernel_run.head(kernel_length).setConstant(1.0 / (2*(double)m));
//   kernel_run(kernel_length-1) *= 0.5;
//   kernel_run.tail(kernel_length-1).setConstant(1.0 / (2*(double)m));
//   kernel_run(n_kernel*2 - kernel_length + 2) *= 0.5;
//
//   if(n_spans > 1) {
//
//     for (size_t i = 1; i < n_spans; ++i) {
//
//       kernel.setZero();
//       m = spans(i);
//       kernel_length = (m + 1);
//
//       kernel.head(kernel_length).setConstant(1.0 / (2*(double)m));
//       kernel(kernel_length-1) *= 0.5;
//       kernel.tail(kernel_length-1).setConstant(1.0 / (2*(double)m));
//       kernel(n_kernel*2 - kernel_length + 2) *= 0.5;
//
//
//       kernel_run = convolve_vec(kernel_run, kernel);
//
//     }
//   }
//
//   out.head(n_kernel + 1) = kernel_run.head(n_kernel + 1);
//   out.tail(n_kernel) = kernel_run.tail(n_kernel);
//
//   return(out);
// }
// //==============================================================================

//==============================================================================
//' @title
//' modified_daniell
//'
//' @description
//' Create a modified daniell kernel using FFT. Adapted from `spec.pgram`.
//'
//' @inheritParams spec.pgram
//'
//' @param n length of the kernel (integer)
//'
//' @return modified Daniell kernel.
//'
//' @noRd
//'
// [[Rcpp::export]]
Eigen::VectorXd modified_daniell(Eigen::VectorXi spans) {

  size_t n_kernel = spans.sum();
  size_t n_spans = spans.size();
  size_t m, kernel_length;


  Eigen::VectorXd kernel(n_kernel * 2 + 1);
  kernel.setZero();
  Eigen::VectorXd kernel_run(n_kernel * 2 + 1);
  kernel_run.setZero();

  m = spans(0);
  kernel_length = (m + 1);

  kernel_run.segment(m, kernel_length).setConstant(1.0 / (2 * (double)m));
  kernel_run(m + kernel_length - 1) *= 0.5;
  kernel_run.head(kernel_length - 1) = kernel_run.segment(m + 1, kernel_length - 1).reverse();

  if(n_spans > 1) {

    for (size_t i = 1; i < n_spans; ++i) {

      kernel.setZero();
      m = spans(i);
      kernel_length = (m + 1);


      kernel.segment(m, kernel_length).setConstant(1.0 / (2 * (double)m));
      kernel(m + kernel_length - 1) *= 0.5;
      kernel.head(kernel_length - 1) = kernel.segment(m + 1, kernel_length - 1).reverse();

      kernel_run = convolve_vec(kernel_run, kernel);

    }
  }

  return(kernel_run);
}

// //==============================================================================
// //' @title
// //' kernel_apply
// //'
// //' @description
// //' Create a modified daniell kernel using FFT. Adapted from `spec.pgram`. This
// //' only calculates the upper triangle when truncated is FALSE.  When truncated
// //' is TRUE the first row is skipped.
// //'
// //' @inheritParams spec.pgram
// //'
// //'
// //' @return modified Daniell kernel.
// //'
// //' @noRd
// //'
// // [[Rcpp::export]]
// Eigen::MatrixXcd kernel_apply(Eigen::MatrixXcd& x,
//                               Eigen::VectorXd& y) {
//
//   Eigen::FFT<double> fft1;
//   Eigen::FFT<double> fft2;
//
//   size_t n_y = y.size();
//   size_t n_row  = x.rows();
//
//   if(n_y != n_row) {
//     Rcpp::stop("kernel_apply: the number of rows of x should equal the length of y");
//   }
//
//   size_t n_col = x.cols();
//   Eigen::VectorXcd x_fft(n_row);
//   Eigen::VectorXcd out_vec(n_row);
//   Eigen::MatrixXcd out(n_row, x.cols());
//
//   Eigen::VectorXcd kernel_fft(n_row);
//   fft1.fwd(kernel_fft, y, 0);
//
//   for (size_t i = 0; i < n_col; ++i) {
//
//     fft2.fwd(x_fft, x.col(i), 0);
//     x_fft = x_fft.array() * kernel_fft.array();
//     fft2.inv(out_vec, x_fft);
//     out.col(i) = out_vec;
//
//   }
//
//
//   return(out);
//
//   // truncated 3x3
//   // -,-,-
//   // 3,4,5
//   // 6,-,8
//
//   // non-truncated 3x3
//   // 0,1,2
//   // -,4,5
//   // -,-,8
//
//   // truncated 4x4
//   //  -, -, -, -
//   //  4, 5, 6, 7
//   //  8, -,10,11
//   // 12, -, -,15
//
//   // non-truncated 4x4
//   //  0, 1, 2, 3
//   //  -, 5, 6, 7
//   //  -, -,10,11
//   //  -, -, -,15
//   //
// }
// //==============================================================================


//==============================================================================
//' @title
//' kernel_apply
//'
//' @description
//' Create a modified daniell kernel using FFT. Adapted from `spec.pgram`. This
//' only calculates the upper triangle when truncated is FALSE.  When truncated
//' is TRUE the first row is skipped.
//'
//' @inheritParams spec.pgram
//'
//'
//' @return modified Daniell kernel.
//'
//' @noRd
//'
// [[Rcpp::export]]
Eigen::MatrixXcd kernel_apply(Eigen::MatrixXcd& x,
                               Eigen::VectorXd& y) {

  size_t n_x = x.rows();
  size_t n_col = x.cols();
  size_t n_y = y.size();
  size_t n_out = n_y * 2 + n_x;
  size_t y_half = n_y / 2;

  VectorXcd wrap(n_y * 2 + n_x);
  MatrixXcd out(n_x, n_col);

  RowVectorXd y_t = y.reverse().transpose();

  for(size_t j=0; j < n_col; ++j) {

    wrap << x.col(j).tail(n_y), x.col(j), x.col(j).head(n_y);

    RcppThread::parallelFor(0, n_x, [&] (size_t i) {

      double re = 0;
      double im = 0;

      re = y_t * wrap.segment(i+2, n_y).real();
      im = y_t * wrap.segment(i+2, n_y).imag();
      out(i, j) = std::complex<double>(re, im);

    });
  }

  return(out);

}
//==============================================================================

//==============================================================================
//' @title
//' spec_taper
//'
//' @description
//' Create a cosine-bell taper. Adapted from `spec.taper`.
//'
//' @inheritParams spec.taper
//'
//'
//' @return cosine-bell taper.
//'
//' @noRd
//'
// [[Rcpp::export]]
Eigen::VectorXd spec_taper(size_t n_row, double p) {


  if(p > 0.5) {
    Rcpp::stop("spec_taper: p cannot be greater than 0.5");
  }
  if(p < 0.0) {
    Rcpp::stop("spec_taper: p cannot be less than 0.0");
  }

  Eigen::VectorXd out(n_row);
  size_t m = floor(n_row * p);

  Eigen::VectorXd taper(n_row);
  taper.setConstant(1.0);

  Eigen::VectorXd ends = M_PI * VectorXd::LinSpaced(m, 1, 2 * m - 1) / (2 * m);

  size_t n_taper = ends.size();
  ends = 0.5 * (1.0 - ends.array().cos() );

  taper.head(n_taper) = ends;
  taper.tail(n_taper) = ends.reverse();

  return(taper);

}
//==============================================================================

//==============================================================================
// might want a better way to get the length and spacings perfect to avoid
// negative values or a smaller than desired last interval size.
//' @title
//' make_groups
//'
//' @description
//' Create a set of lengths with increasing sizes useful from grouping results.
//' The first and last groups have length equal to 1.
//'
//' @param n the number of values to subset into groups. (integer)
//' @param max_lag how fast groups get bigger. Larger numbers have a larger range
//' of group sizes. (integer)
//' @param n_groups the number of groups to create. (integer)
//' @param min_aggregate the minimum size for a group. (integer)
//'
//'
//' @return Matrix with filled in complex conjugate columns.
//'
//' @noRd
//'
// [[Rcpp::export]]
Eigen::VectorXi make_groups(size_t n_groups, size_t n) {


  // check inputs
  if (n_groups <= 0) {
    Rcpp::stop("log_lags_eigen: n must be greater than 0");
  }


  if (n < 0) {
    Rcpp::stop("log_lags_eigen: max_time_lag must be non-negative");
  }


  if(n_groups > (n + 1L)) {
    Rcpp::warning("The number of lags is greater than the maximum time lag");
    return(VectorXi::LinSpaced(n + 1, 0, n));
  }


  // Lags begin at zero
  ArrayXd lags = pow(10.0, ArrayXd::LinSpaced(n_groups, 0.0, std::log10(n))) - 1.0;
  lags = lags.round();


  // difference between values
  lags.head(n_groups - 1) = lags.tail(n_groups - 1) - lags.head(n_groups - 1);
  lags = lags.array().max(1.0);
  std::sort(lags.data(), lags.data() + lags.size());


  lags(0) = 1;
  lags(lags.size() - 1) = 1;
  size_t s = lags.sum() - n;


  if(s > 0) {
    for(size_t i = lags.size() - 2; i > lags.size() - 2 - s; --i) {
      lags(i) = lags(i) - 1.0;
    }
  }


  VectorXi lg(lags.size());
  for(size_t i = 0; i < lags.size(); ++i) {
    lg(i) = int(lags(i));
  }

  return(lg);
}

// //==============================================================================
// // might want a better way to get the length and spacings perfect to avoid
// // negative values or a smaller than desired last interval size.
// //' @title
// //' make_groups
// //'
// //' @description
// //' Create a set of lengths with increasing sizes useful from grouping results.
// //' The first and last groups have length equal to 1.
// //'
// //' @param n_row the number of values to subset into groups. (integer)
// //' @param power how fast groups get bigger. Larger numbers have a larger range
// //' of group sizes. (integer)
// //' @param n_groups the number of groups to create. (integer)
// //' @param min_aggregate the minimum size for a group. (integer)
// //'
// //'
// //' @return Matrix with filled in complex conjugate columns.
// //'
// //' @noRd
// //'
// // [[Rcpp::export]]
// Eigen::VectorXi make_groups(size_t n_row,
//                             size_t power,
//                             size_t n_groups,
//                             size_t min_aggregate) {
//
//   if (n_groups >= n_row) {
//     Rcpp::warning("make_groups: min_aggregate is ignored when n_groups is greater than n_row. Returning a group for each row.");
//     return(Eigen::VectorXi::Ones(n_row));
//   }
//
//   if ( ((n_row - 2) / n_groups) < min_aggregate) {
//     Rcpp::stop("make_groups: min_aggregate is too large for given n_row and n_groups");
//   }
//
//
//   Eigen::VectorXd groups_f = VectorXd::LinSpaced(n_groups, 0, 1);
//   groups_f = groups_f.array().pow(power);
//   groups_f(0) = 0.0;
//   groups_f(n_groups - 1) = 0.0;
//   double scale = std::floor(n_row / groups_f.sum());
//
//   groups_f = ((groups_f * scale).array().max((double)min_aggregate).round());
//   groups_f(0) = 1.0;
//   groups_f(n_groups - 1) = 1.0;
//
//   // create the group sizes
//   // first and last group are DC components not to blend
//   Eigen::VectorXi groups(n_groups);
//   for (size_t i = 0; i < n_groups; ++i) {
//     groups(i) = (size_t)groups_f(i);
//   }
//
//   size_t g_sum = groups.sum();
//
//   // decrease last non-DC group if sum is too large
//   if(g_sum > n_row) {
//
//     // updating/improving the method could decrease the number of times this
//     // statement is entered
//     if((g_sum - n_row) >=  groups(n_groups - 2)) {
//       Rcpp::stop("make_groups: The number of groups, power or min_aggregate is too large. Try reducing power, min_aggregate or the number of groups.");
//     }
//
//     groups(n_groups-2) -= (g_sum - n_row);
//   }
//
//   // increase last non-DC group if sum is too small
//   if(g_sum < n_row) {
//     groups(n_groups-2) += (n_row - g_sum);
//   }
//
//   return(groups);
//
// }
// //==============================================================================


//==============================================================================
//' @title
//' power_spaced
//'
//' @description
//' Create an n length vector of power spaced series between min and max.
//'
//' @param n the number of values. (integer)
//' @param min minimum value in series. (integer)
//' @param max maximum value in series. (integer)
//' @param power how fast values change. (integer)
//'
//'
//' @return vector of power spaced series between min and max.
//'
//' @noRd
//'
// [[Rcpp::export]]
Eigen::ArrayXd power_spaced(size_t n, double min, double max, double power) {

  Eigen::ArrayXd out = ArrayXd::LinSpaced(n, 0, 1);

  out = out.pow(power);

  out = out * (max - min) + min;

  return(out);

}
//==============================================================================


//==============================================================================
//' @title
//' group_frequency
//'
//' @description
//' Get the mean frequency for each group.
//'
//' @inheritParams make_groups
//'
//' @param frequency vector of frequencies. (integer)
//'
//'
//' @return vector frequencies for each group.
//'
//' @noRd
//'
// [[Rcpp::export]]
Eigen::VectorXd group_frequency(Eigen::ArrayXd frequencies,
                                size_t n_groups) {

  Eigen::VectorXd fs(n_groups);
  size_t s = 0; // start of the group
  size_t group_size = 0;

  Eigen::VectorXi g = make_groups(n_groups, frequencies.size());

  for(size_t i = 0; i < n_groups; ++i) {
    group_size = g(i);

    fs(i) = frequencies.segment(s, group_size).array().mean();
    s += group_size;

  }

  return(fs);

}
//==============================================================================


//==============================================================================
//' @title
//' determine_frequency
//'
//' @description
//' Get frequencies based on series length
//'
//'
//' @param n length of series. (integer)
//'
//'
//' @return vector of frequencies between 0 and 0.5.
//'
//' @noRd
//'
// [[Rcpp::export]]
Eigen::VectorXd determine_frequency(size_t n) {

  return(Eigen::VectorXd::LinSpaced((double)(n / 2), 0, 0.5 - 1.0 / (double)n));

}
//==============================================================================


//==============================================================================
// currently unused
// [[Rcpp::export]]
Eigen::MatrixXcd check_ffts(Eigen::MatrixXcd& x,
                            double cutoff) {

  std::complex<double> cplx_zero = std::complex<double>(1e12, 1e12);

  double ratio = 0;

  size_t n_row = x.rows();
  size_t n_col = x.cols();

  Eigen::MatrixXd x_abs = x.cwiseAbs2();
  Eigen::VectorXd x_max(n_col);

  for(size_t j = 0; j < n_col; ++j) {
    x_max(j) = x_abs.col(j).maxCoeff();
    x_abs.col(j) = x_abs.col(j).array() / x_max(j);
  }

  for(size_t i = 0; i < n_row; ++i) {
    for(size_t j = 1; j < n_col; ++j) {

      ratio = x_abs(i, 0) / x_abs(i, j);
      // ratio = x_max(j) / x_abs(i, j);

      if (ratio > cutoff) {
        x(i, j) = cplx_zero;
      }

    }

  }

  return(x);
}
//==============================================================================



//==============================================================================
//' @title
//' which_indices
//'
//' @description
//' Determine the intervals that each x value falls in
//'
//' @inheritParams splines::bs
//' @param knots location of knots for the b-splines. Unlike `splines::bs` this
//' includes the boundary knots. (numeric vector)
//'
//'
//' @return the interval ids that each x value falls in.
//'
//' @noRd
//'
// [[Rcpp::export]]
Eigen::VectorXi which_indices(const Eigen::VectorXd& x,
                              const Eigen::VectorXd& knots) {

  size_t n_x = x.size();
  size_t n_below = 0;
  size_t n_old = 0;
  size_t n_i = 0;
  size_t n = 0;
  size_t n_knots = knots.size();
  Eigen::VectorXi ind = Eigen::VectorXi::Constant(n_x, n_knots - 2);


  for (size_t i = 0; i < n_knots; ++i) {

    n_i = (x.array() <= knots(i)).count();
    n_below = n_i - n_old;

    if(n_below > 0) {
      if (i == 0) {
        ind.segment(n, n_below) = Eigen::VectorXi::Constant(n_below, i);
      } else {
        ind.segment(n, n_below) = Eigen::VectorXi::Constant(n_below, i - 1);
      }
    }

    n_old = n_i;
    n += n_below;

  }

  return(ind);
}
//==============================================================================


/*** R
*/
