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

  VectorXd out;

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


  VectorXd hann = window_hann(n);
  VectorXcd out(n);

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



// =============================================================================
//' @title
//' window_first_deriv
//'
//' @description
//' First derivative window for FFT
//'
//' @inheritParams window_hann
//'
//' @param a0 \code{double} coefficient
//' @param a1 \code{double} coefficient
//' @param a2 \code{double} coefficient
//' @param a3 \code{double} coefficient
//'
//' @return window
//'
//' @export
//'
//' @examples
//' # nuttall window
//' window_first_deriv(100, 0.355768, 0.487396, 0.144232, 0.012604)
//'
//' @noRd
//'
// [[Rcpp::export]]
Eigen::ArrayXd window_first_deriv(size_t n,
                                   double a0,
                                   double a1,
                                   double a2,
                                   double a3) {

  ArrayXd k = Eigen::ArrayXd::LinSpaced(n, 0.0, n - 1.0) / double(n - 1.0);

  return(a0 - a1 * (2 * M_PI * k).cos() +
              a2 * (4 * M_PI * k).cos() -
              a3 * (6 * M_PI * k).cos());

}

// =============================================================================
//' @title
//' window_nuttall
//'
//' @description
//' Nuttall window for FFT
//'
//' @inheritParams window_hann
//'
//' @return window
//'
//' @export
//'
//' @examples
//' window_nuttall(100)
//'
//' @noRd
//'
// [[Rcpp::export]]
Eigen::ArrayXd window_nuttall(size_t n) {

  double a0 = 0.355768;
  double a1 = 0.487396;
  double a2 = 0.144232;
  double a3 = 0.012604;

  return(window_first_deriv(n, a0, a1, a2, a3));

}


// =============================================================================
//' @title
//' window_blackman_nuttall
//'
//' @description
//' Blackman-Nuttall window for FFT
//'
//' @inheritParams window_hann
//'
//' @return window
//'
//' @export
//'
//' @examples
//' window_blackman_nuttall(100)
//'
//' @noRd
//'
// [[Rcpp::export]]
Eigen::ArrayXd window_blackman_nuttall(size_t n) {

  double a0 = 0.3635819;
  double a1 = 0.4891775;
  double a2 = 0.1365995;
  double a3 = 0.0106411;

  return(window_first_deriv(n, a0, a1, a2, a3));

}

// =============================================================================
//' @title
//' window_blackman_harris
//'
//' @description
//' Blackman-Harris window for FFT
//'
//' @inheritParams window_hann
//'
//' @return window
//'
//' @export
//'
//' @examples
//' window_blackman_harris(100)
//'
//' @noRd
//'
// [[Rcpp::export]]
Eigen::ArrayXd window_blackman_harris(size_t n) {

  double a0 = 0.35875;
  double a1 = 0.48829;
  double a2 = 0.14128;
  double a3 = 0.01168;

  return(window_first_deriv(n, a0, a1, a2, a3));

}

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
