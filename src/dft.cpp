#include "hydrorecipes.h"


// [[Rcpp::export]]
std::complex<double> dft(Eigen::VectorXd x,
                         double frequency) {

  unsigned int n = x.size();
  double n_2 = (double)n / 2;

  Eigen::VectorXd k = Eigen::VectorXd::LinSpaced(n, 0.0, (double)n - 1.0);
  std::complex<double> two_i(0.0, 2.0);
  std::complex<double> dft =
    (x.array() * (-two_i * M_PI * k.array() * frequency / n).exp()).sum() / n_2;

  return(dft);
}

// [[Rcpp::export]]
std::complex<double> dft_with_window(Eigen::VectorXd x,
                                     double frequency) {
  unsigned int n = x.size();
  Eigen::VectorXd window = window_hann(n);

  x = x.array() * window.array();

  double n_2 = (double)n / 2;

  Eigen::VectorXd k = Eigen::VectorXd::LinSpaced(n, 0.0, (double)n - 1.0);
  std::complex<double> two_i(0.0, 2.0);
  std::complex<double> dft =
    (x.array() * (-two_i * M_PI * k.array() * frequency / n).exp()).sum() / n_2;

  return(dft);
}

// [[Rcpp::export]]
std::complex<double> dft_goertzel(Eigen::VectorXd x,
                                  double frequency) {


  unsigned int n = x.size();
  double scale = 0.5 * (double)n;
  Eigen::VectorXd window = window_hann(n);

  x = x.array() * window.array();

  double omega = 2.0 * M_PI * frequency / (double)n;

  double sin_omega = sin(omega); // imageW
  double cos_omega = cos(omega);

  double co = 2.0 * cos_omega; // realW
  double term_0 = 0.0;
  double term_1 = 0.0;
  double term_2 = 0.0;

  for (unsigned int i = 0; i < n; ++i) {
    term_0 = x[i] + co * term_1 - term_2;
    term_2 = term_1;
    term_1 = term_0;
  }

  return(std::complex<double>((term_1 - cos_omega * term_2) / scale,
                              (term_1 * sin_omega) / scale));
}


// [[Rcpp::export]]
Eigen::MatrixXcd be_dft(Eigen::MatrixXd x,
                        double frequency) {

  unsigned int n = x.rows();
  double n_2 = (double)n / 2;

  Eigen::MatrixXcd fft_matrix(1,3);

  x.col(0) = detrend_and_demean_matrix(x.col(0), true, true);
  x.col(1) = detrend_and_demean_matrix(x.col(1), true, true);
  x.col(2) = detrend_and_demean_matrix(x.col(2), true, true);


  fft_matrix(0) = dft_with_window(x.col(0), frequency);
  fft_matrix(1) = dft_with_window(x.col(1), frequency);
  fft_matrix(2) = dft_with_window(x.col(2), frequency);

  Eigen::MatrixXcd dft_mat = multiply_ffts(fft_matrix);


  return(solve_cplx_parallel(dft_mat));

}
