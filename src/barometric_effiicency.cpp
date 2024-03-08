#include "frecipes.h"


//' @title
//' be_clark_cpp
//'
//' @description
//' Clark 1967 solution for calculating barometric efficiency (Algorithm from Batu 1998, pg 76)
//'
//' @param dep \code{numeric vector} of the dependent variable (ie:water level)
//' @param ind \code{numeric vector} of the independent variable (ie:barometric pressure)
//' @param lag_space \code{integer} spacing for lags, useful for higher frequency monitoring
//' @param inverse \code{logical} whether the barometric relationship is inverse
//'
//' @return barometric efficiency using Clark's method
//'
//'
//' @export
//'
//' @examples
//' n <- 1000
//' baro <- sin(seq(0, 2*pi, length.out = 1000))
//' wl <- -0.4 * baro + rnorm(1000, sd = 0.02)
//' be_clark_cpp(wl, baro, lag_space=1, inverse=TRUE)
//'
// [[Rcpp::export]]
double be_clark_cpp(arma::vec dep,
                    arma::vec ind,
                    int lag_space,
                    bool inverse) {

  int n = dep.n_elem - lag_space;
  arma::vec ret;


  dep = dep.head( n ) - dep.tail( n );
  ind = ind.head( n ) - ind.tail( n );

  if (inverse) {
    dep = -dep;
  }

  dep = cumsum((arma::sign(dep) % arma::sign(ind)) % abs(dep));
  ind = arma::cumsum(abs(ind));

  ret = arma::solve(ind, dep);

  return ret(0);
}


// // [[Rcpp::export]]
// double be_clark_eigen(Eigen::MatrixXd y,
//                       Eigen::MatrixXd x,
//                       unsigned int lag_space,
//                       bool inverse) {
//
//   unsigned int n_y = y.size();
//   unsigned int n_x = x.size();
//   unsigned int n = n_y - lag_space;
//
//   if (n_x != n_y) {
//     Rcpp::stop("the lengths of x and y must be equal");
//   }
//
//   if (n_y <= n) {
//     Rcpp::stop("lag_space cannot be equal to or larger than the length of x and y");
//   }
//
//
//   Eigen::MatrixXd dep = y.head(n).array() - y.tail(n).array();
//   Eigen::MatrixXd ind = ind.head(n).array() - x.tail(n).array();
//
//   if (inverse) {
//     dep.array() = -dep.array();
//   }
//
//   dep = ((dep.array().sign() * ind.array().sign()) * dep.array().abs());
//   ind = (ind.array().abs());
//
//   for (unsigned int i = 0; i < ind.size() - 1; ++i) {
//     dep[i+1] = dep[i] + dep[i+1];
//     ind[i+1] = ind[i] + ind[i+1];
//   }
//
//   Eigen::MatrixXd ind_m = Eigen::MatrixXd::Map(ind.data(), n, 1);
//   Eigen::MatrixXd dep_m = Eigen::MatrixXd::Map(dep.data(), n, 1);
//
//   Eigen::MatrixXd out = llt_solve(ind_m, dep_m);
//
//   // Eigen::LLT<Eigen::MatrixXd> llt(Eigen::MatrixXd(ind).setZero().selfadjointView<Lower>().
//   //                                      rankUpdate(ind.adjoint()));
//   // Eigen::MatrixXd out = (llt.solve(ind.adjoint() * dep));
//
//
//   // Eigen::VectorXd out = ind.colPivHouseholderQr().solve(dep);
//   return(out[0]);
// }


// [[Rcpp::export]]
double be_least_squares_diff_cpp(arma::vec dep,
                            arma::vec ind,
                            int lag_space,
                            bool inverse) {

  int n = dep.n_elem-lag_space;
  arma::vec ret;


  dep = dep.head( n ) - dep.tail( n );
  ind = ind.head( n ) - ind.tail( n );

  if (inverse) {
    dep = -dep;
  }

  ret = arma::solve(ind, dep);

  return ret(0);
}


// [[Rcpp::export]]
double be_least_squares_cpp(arma::vec dep,
                            arma::vec ind,
                            bool inverse) {

  unsigned int n = ind.n_elem;
  arma::mat y = arma::ones(n, 2);
  y.col(0) = ind;

  arma::vec ret;

  if (inverse) {
    dep = -dep;
  }

  ret = arma::solve(y, dep);

  return ret(0);
}

//' @title
//' be_acworth_calc_cpp
//'
//' @description
//' Acworth, R. I., Halloran, L. J. S., Rau, G. C., Cuthbert, M. O.,
//'  & Bernardi, T. L. (2016). An objective frequency-domain method for
//'   quantifying confined aquifer compressible storage using Earth and
//'   atmospheric tides. Geophysical Research Letters, 43(November).
//'   https://doi.org/10.1002/2016GL071328
//'
//' @param s2_gw \code{numeric} s2 component in the groundwater levels
//' @param s2_et \code{numeric} s2 component in the earth tides
//' @param s2_at \code{numeric} s2 component for atmospheric pressure
//' @param m2_gw \code{numeric} m2 component in the groundwater levels
//' @param m2_et \code{numeric} m2 component in the earth tides
//' @param d_phase \code{numeric} phase difference between Earth tide and atmospheric drivers s2_et and s2_at
//' @param inverse \code{logical} whether the barometric relationship is inverse (TRUE means that when the barometric pressure goes up the measured water level goes down (vented transducer, depth to water), FALSE means that when the barometric pressure goes up so does the measured pressure (non-vented transducer))
//'
//' @return barometric efficiency
//' @export
//'
//' @examples
//'
//' be_acworth_cpp(s2_at = 7.461,
//'            s2_et=224.640,
//'            s2_gw=4.086,
//'            m2_gw = 0.471,
//'            m2_et = 492.526,
//'            d_phase=-56.709,
//'            inverse = TRUE)
//' be_acworth_cpp(s2_at = 6.164,
//'            s2_et=270.463,
//'            s2_gw=0.329,
//'            m2_gw = 0.225,
//'            m2_et = 551.572,
//'            d_phase=-71.726,
//'            inverse = TRUE)
//' be_acworth_cpp(s2_at = 5.897,
//'            s2_et=234.478,
//'            s2_gw=5.536,
//'            m2_gw = 0.773,
//'            m2_et = 558.075,
//'            d_phase=-70.393,
//'            inverse = TRUE)
// [[Rcpp::export]]
double be_acworth_calc_cpp(const double s2_gw,
                           const double s2_et,
                           const double s2_at,
                           const double m2_gw,
                           const double m2_et,
                           const double d_phase,
                           const bool inverse) {


  double term = 0.0;

  // check if earthtides are present
  if ((m2_et != 0) & (s2_et != 0)) {
    term = s2_et * (cos(d_phase) * m2_gw / m2_et);
  }

  Rcpp::Rcout << "term: " << term << std::endl;

  // inverse equal to TRUE suggests barometric efficiency
  if (inverse) {
    return ((s2_gw + term) / s2_at);
  }

  // inverse equal to FALSE suggests loading efficiency therefore need to
  // correct
  return (1.0 - (s2_gw - term) / s2_at);
}



// [[Rcpp::export]]
Eigen::Vector2i get_peaks(Eigen::VectorXd freqs, double f1, double f2) {

  if (f1 < freqs.minCoeff()) {
    Rcpp::stop("Chosen frequency (f1) is outside fft range of frequencies");
  }
  if (f2 < freqs.minCoeff()) {
    Rcpp::stop("Chosen frequency (f2) is outside fft range of frequencies");
  }

  unsigned int n = freqs.size();
  unsigned int done = 0;
  Eigen::Vector2i r = Eigen::Vector2i::Zero();

  for (unsigned int j = 0; j < (n-1); ++j) {
    if (done == 2) {
      break;
    }

    if ((f1 >= freqs[j]) & (f1 <= freqs[j + 1])) {
      done += 1;
      if (std::abs(f1 - freqs[j]) < std::abs(f1 - freqs[j + 1])) {
        r[0] = j;
      } else {
        r[0] = j + 1;
      }
    }

    if ((f2 >= freqs[j]) & (f2 <= freqs[j + 1])) {
      done += 1;
      if (std::abs(f2 - freqs[j]) < std::abs(f2 - freqs[j + 1])) {
        r[1] = j;
      } else {
        r[1] = j + 1;
      }
    }

  }

  return(r);
}


// // [[Rcpp::export]]
// double complex_mod(std::complex<double> x) {
//   double x_r = (double)x.real();
//   double x_i = (double)x.imag();
//   return sqrt(x_r * x_r + x_i * x_i);
// }

// [[Rcpp::export]]
double be_acworth_cpp(Eigen::MatrixXd& x,
                      const Eigen::VectorXi& spans,
                      bool detrend,
                      bool demean,
                      double taper,
                      bool inverse,
                      double f1,
                      double f2,
                      double frequency_scale) {



  Eigen::MatrixXcd pgram = spec_pgram(x,
                                      spans,
                                      detrend,
                                      demean,
                                      taper);

  unsigned int n = x.rows();
  Eigen::VectorXd freqs = determine_frequency(n) * frequency_scale;
  Eigen::Vector2i r = get_peaks(freqs, f1, f2);

  // get order of frequencies
  unsigned int s = 1;
  unsigned int m = 0;
  if (r[0] > r[1]) {
    s = 0;
    m = 1;
  }

  double s2_gw, m2_gw, s2_at, m2_at, s2_et, m2_et, d_phase;

  s2_gw = std::sqrt(pgram(r[s], 0).real());
  m2_gw = std::sqrt(pgram(r[m], 0).real());

  s2_at = std::sqrt(pgram(r[s], 3).real());
  m2_at = std::sqrt(pgram(r[m], 3).real());

  s2_et = std::sqrt(pgram(r[s], 5).real());
  m2_et = std::sqrt(pgram(r[m], 5).real());

  d_phase = std::arg(pgram(r[s], 5)) - std::arg(pgram(r[s], 3));


  double be = be_acworth_calc_cpp(s2_gw,
                                  s2_et,
                                  s2_at,
                                  m2_gw,
                                  m2_et,
                                  d_phase,
                                  inverse);

  Rcpp::Rcout << "rows: " << pgram.rows() << std::endl;
  Rcpp::Rcout << "cols: " << pgram.cols() << std::endl;

  Rcpp::Rcout << "s2_gw: " << s2_gw << std::endl;
  Rcpp::Rcout << "m2_gw: " << m2_gw << std::endl;
  Rcpp::Rcout << "s2_at: " << s2_at << std::endl;
  Rcpp::Rcout << "m2_at: " << m2_at << std::endl;
  Rcpp::Rcout << "s2_et: " << s2_et << std::endl;
  Rcpp::Rcout << "m2_et: " << m2_et << std::endl;

  return(be);

}



