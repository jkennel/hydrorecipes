#include "hydrorecipes.h"

//==============================================================================
//' @title
//' window_hann
//'
//' @description
//' Hann window for FFT.
//'
//' @param n length of the window vector (integer)
//'
//' @return window of length n.
//'
//' @noRd
//'
// [[Rcpp::export]]
Eigen::VectorXd window_hann(size_t n) {

  Eigen::VectorXd out;

  if(n == 0) {
    Rcpp::stop("window_hann: n must be an integer greater than 0");
  }
  if (n == 1) {
    out = VectorXd::Ones(1);
    return(out);
  }
  if (n == 2) {
    out = VectorXd::Ones(2) / 2.0;
    return(out);
  }

  out = Eigen::VectorXd::LinSpaced(n, 0.0, (double)n - 1.0);
  double denom = 1.0 / (n - 1);
  out = M_2PI * out.array() * denom;
  out = 0.5 * (1.0 - out.array().cos());

  return(out);

}
//==============================================================================


//==============================================================================
//' @title
//' window_tukey
//'
//' @description
//' Tukey window for FFT.
//'
//' @inheritParams window_hann
//' @param r percent on each side to taper (double)
//'
//' @return window of length n.
//'
//' @noRd
//'
// [[Rcpp::export]]
Eigen::VectorXd window_tukey(size_t n, double r) {

  size_t n_lobe = std::ceil(r * (double)n / 2.0);

  if(n == 0) {
    Rcpp::stop("window_tukey: n must be an integer greater than 0");
  }
  if (n == 1) {
    return(VectorXd::Ones(1));
  }
  if (n == 2) {
    return(VectorXd::Zero(2));
  }

  if(r <= 0) {
    return(VectorXd::Ones(n));
  } else if (r >= 1) {
    return(window_hann(n));
  }

  double mult = 1.0 / ((double)n - 1.0) * ((double)n_lobe - 1.0) * M_2PI  / r;

  VectorXd i_vals = VectorXd::LinSpaced(n_lobe, 0.0, mult);
  VectorXd out = VectorXd::Ones(n);
  VectorXd ends = 0.5 * (1.0 - i_vals.array().cos());

  out.head(n_lobe) = ends;
  out.tail(n_lobe) = ends.tail(n_lobe).reverse();


  return(out);
}

//==============================================================================
//' @title
//' window_hann_cplx
//'
//' @description
//' Hann window for complex FFT.
//'
//' @inheritParams window_hann
//'
//' @return window of length n.
//'
//' @noRd
//'
// [[Rcpp::export]]
Eigen::VectorXcd window_hann_cplx(size_t n) {


  Eigen::VectorXd hann = window_hann(n);
  Eigen::VectorXcd out(n);

  for(size_t i = 0; i < out.size(); ++i) {
    out(i) = std::complex<double>(hann(i), hann(i));
  }

  return(out);

}
//==============================================================================


//==============================================================================
//' @title
//' window_rectangle
//'
//' @description
//' Rectangular window.
//'
//' @inheritParams window_hann
//'
//' @return window of length n.
//'
//' @noRd
//'
// [[Rcpp::export]]
Eigen::VectorXd window_rectangle(size_t n) {

  return(VectorXd::Ones(n)/(double)n);

}
//==============================================================================


//==============================================================================
//' @title
//' window_scale
//'
//' @description
//' Scale factor for a window function.
//'
//' @param window the window function (numeric vector)
//' @param n_new length of the padded series (integer)
//' @param n_fft length of the input series (integer)
//'
//' @return window of length n.
//'
//' @noRd
//'
// [[Rcpp::export]]
double window_scale(Eigen::VectorXd window,
                    size_t n_new,
                    size_t n_fft) {

  double window_sum  = window.array().sum();
  double window_mean = window.array().mean();
  double window_mean_sq = window_mean * window_mean;
  double window_sum_sq = window_sum * window_sum;
  // double window_sq_sum = (window.array() * window.array()).sum();
  double window_sq_sum = (window.array().square()).sum();
  double scale = window_mean_sq * ((double)n_new * window_sq_sum / window_sum_sq) * (double)n_new * (double)n_fft;

  scale = 1.0 / scale;

  return(scale);
}
//==============================================================================



/*** R
*/
