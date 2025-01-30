#include "hydrorecipes.h"


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
//' @export
//'
//' @examples
//' n <- 1000
//' baro <- sin(seq(0, 2 * pi, length.out = 1000))
//' wl <- -0.4 * baro + rnorm(1000, sd = 0.02)
//' be_clark_cpp(wl, baro, lag_space = 1, inverse = TRUE)
//'
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

  dep = arma::cumsum((arma::sign(dep) % arma::sign(ind)) % abs(dep));
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
Eigen::MatrixXd be_least_squares_diff_cpp(Eigen::VectorXd dep,
                                 Eigen::VectorXd ind,
                                 int lag_space,
                                 bool inverse) {

  unsigned int n = dep.size() - lag_space;

  Eigen::MatrixXd x = Eigen::MatrixXd::Ones(n, 1);
  Eigen::MatrixXd y = Eigen::MatrixXd(n, 1);

  if (inverse) {
    dep = -dep;
  }

  // difference inputs
  y.col(0) = dep.head( n ) - dep.tail( n );
  x.col(0) = ind.head( n ) - ind.tail( n );

  const int p = 1;

  const Eigen::LLT<Eigen::MatrixXd> llt(Eigen::MatrixXd(p, p).setZero().selfadjointView<Lower>().
                                          rankUpdate(x.adjoint()));

  return llt.solve(x.adjoint() * y);

}


// [[Rcpp::export]]
Eigen::MatrixXd be_least_squares_cpp(Eigen::VectorXd dep,
                            Eigen::VectorXd ind,
                            bool inverse) {


  unsigned int n = ind.size();
  Eigen::MatrixXd x = Eigen::MatrixXd(n, 1);
  Eigen::MatrixXd y = Eigen::MatrixXd(n, 1);

  if (inverse) {
    dep = -dep;
  }

  // detrend inputs
  y.col(0) = detrend_vector(dep);
  x.col(0) = detrend_vector(ind);

  const int p = 1;

  const Eigen::LLT<Eigen::MatrixXd> llt(Eigen::MatrixXd(p, p).setZero().selfadjointView<Lower>().
                                          rankUpdate(x.adjoint()));


  return llt.solve(x.adjoint() * y);
}

// //' @title
// //' be_acworth_calc_cpp
// //'
// //' @description
// //' Acworth, R. I., Halloran, L. J. S., Rau, G. C., Cuthbert, M. O.,
// //'  & Bernardi, T. L. (2016). An objective frequency-domain method for
// //'   quantifying confined aquifer compressible storage using Earth and
// //'   atmospheric tides. Geophysical Research Letters, 43(November).
// //'   https://doi.org/10.1002/2016GL071328
// //'
// //' @param s2_gw \code{numeric} s2 component in the groundwater levels
// //' @param s2_et \code{numeric} s2 component in the earth tides
// //' @param s2_at \code{numeric} s2 component for atmospheric pressure
// //' @param m2_gw \code{numeric} m2 component in the groundwater levels
// //' @param m2_et \code{numeric} m2 component in the earth tides
// //' @param d_phase \code{numeric} phase difference between Earth tide and atmospheric drivers s2_et and s2_at
// //' @param inverse \code{logical} whether the barometric relationship is inverse (TRUE means that when the barometric pressure goes up the measured water level goes down (vented transducer, depth to water), FALSE means that when the barometric pressure goes up so does the measured pressure (non-vented transducer))
// //'
// //' @return barometric efficiency
// //' @export
// //'
// //' @examples
// //'
// //' be_acworth_calc_cpp(s2_at = 7.461,
// //'            s2_et=224.640,
// //'            s2_gw=4.086,
// //'            m2_gw = 0.471,
// //'            m2_et = 492.526,
// //'            d_phase=-56.709)
// //' be_acworth_calc_cpp(s2_at = 6.164,
// //'            s2_et=270.463,
// //'            s2_gw=0.329,
// //'            m2_gw = 0.225,
// //'            m2_et = 551.572,
// //'            d_phase=-71.726)
// //' be_acworth_calc_cpp(s2_at = 5.897,
// //'            s2_et=234.478,
// //'            s2_gw=5.536,
// //'            m2_gw = 0.773,
// //'            m2_et = 558.075,
// //'            d_phase=-70.393)
// // [[Rcpp::export]]
// double be_acworth_calc_cpp(const std::complex<double> s2_gw,
//                            const std::complex<double> s2_et,
//                            const std::complex<double> s2_at,
//                            const std::complex<double> m2_gw,
//                            const std::complex<double> m2_et,
//                            double d_phase) {
//
//
//   double term = 0.0;
//
//   // check if earthtides are present
//   if ((std::abs(m2_et) != 0) & (std::abs(s2_et) != 0)) {
//     term = std::abs(s2_et) * (std::cos(d_phase) * std::abs(m2_gw) / std::abs(m2_et));
//   }
//
//   return ((std::abs(s2_gw) - term) / std::abs(s2_at));
//
// }
//
//
// //' @title
// //' be_rau_calc_cpp
// //'
// //' @description
// //' Rau, G.C., Cuthbert, M.O., Acworth, R.I. and Blum, P., 2020.
// //' Disentangling the groundwater response to Earth and atmospheric tides
// //' to improve subsurface characterisation. Hydrology and earth system
// //' sciences, 24(12), pp.6033-6046.
// //'
// //' @param s2_gw \code{numeric} s2 component in the groundwater levels
// //' @param s2_et \code{numeric} s2 component in the earth tides
// //' @param s2_at \code{numeric} s2 component for atmospheric pressure
// //' @param m2_gw \code{numeric} m2 component in the groundwater levels
// //' @param m2_et \code{numeric} m2 component in the earth tides
// //' @param amp_ratio \code{numeric} amplitude ratio to account for damping
// //'
// //' @return barometric efficiency
// //' @export
// //'
// //' @examples
// //'
// // [[Rcpp::export]]
// double be_rau_calc_cpp(const std::complex<double> s2_gw,
//                        const std::complex<double> s2_et,
//                        const std::complex<double> s2_at,
//                        const std::complex<double> m2_gw,
//                        const std::complex<double> m2_et,
//                        const double amp_ratio) {
//
//
//   double term = 0.0;
//
//   // equation 9
//   return(1.0 / amp_ratio * std::abs((s2_gw - (m2_gw / m2_et) * s2_et) / s2_at));
//
//
// }



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

  for (unsigned int j = 0; j < (n - 1); ++j) {
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
// double be_acworth_cpp(Eigen::MatrixXd& x,
//                       const Eigen::VectorXi& spans,
//                       bool detrend,
//                       bool demean,
//                       double taper,
//                       double f1,
//                       double f2,
//                       double frequency_scale) {
//
//
//
//   if (spans.size() < 1) {
//     Rcpp::stop("spec_pgram: spans must be length 1 or larger.");
//   }
//
//   // detrend or demean
//   x = detrend_and_demean_matrix(x, detrend, demean);
//
//   size_t n_row = x.rows();
//   size_t n_new = next_n_eigen(n_row);
//   Rcpp::Rcout << "n_new: " << n_new << std::endl;
//
//   std::complex<double> scale = 1.0 / n_row; // or n_new
//
//   // taper vector
//   ArrayXd taper_array = spec_taper(n_row, taper).array();
//
//
//   // Do FFTs
//   MatrixXcd pgram = fft_matrix(x.array().colwise() * taper_array,
//                                n_row);
//
//
//   unsigned int n = pgram.rows();
//   Eigen::VectorXd freqs = determine_frequency(n) * frequency_scale;
//   Eigen::Vector2i r = get_peaks(freqs, f1, f2);
//
//   Rcpp::Rcout << "r: " << r << std::endl;
//
//   // get order of frequencies
//   unsigned int s = 1;
//   unsigned int m = 0;
//   if (r[0] > r[1]) {
//     s = 0;
//     m = 1;
//   }
//
//   std::complex<double> s2_gw, m2_gw, s2_at, m2_at, s2_et, m2_et;
//
//   s2_gw = pgram(r[s], 0);
//   m2_gw = pgram(r[m], 0);
//
//   s2_at = pgram(r[s], 1);
//   m2_at = pgram(r[m], 1);
//
//   s2_et = pgram(r[s], 2);
//   m2_et = pgram(r[m], 2);
//
//
//   double d_phase = std::arg(s2_at) - std::arg(s2_et);
//
//
//   double be = be_acworth_calc_cpp(s2_gw,
//                                   s2_et,
//                                   s2_at,
//                                   m2_gw,
//                                   m2_et,
//                                   d_phase);
//
//
//   // Rcpp::Rcout << "s2_at: " << std::arg(s2_at) << std::endl;
//   // Rcpp::Rcout << "s2_et: " << std::arg(s2_et) << std::endl;
//   // Rcpp::Rcout << "s2_gw: " << std::arg(s2_gw) << std::endl;
//   // Rcpp::Rcout << "m2_gw: " << std::arg(m2_gw) << std::endl;
//   //
//   // Rcpp::Rcout << "rows: " << pgram.rows() << std::endl;
//   // Rcpp::Rcout << "cols: " << pgram.cols() << std::endl;
//   //
//   // Rcpp::Rcout << "s2_gw: " << 2*std::abs(s2_gw)/n << std::endl;
//   // Rcpp::Rcout << "m2_gw: " << 2*std::abs(m2_gw)/n << std::endl;
//   // Rcpp::Rcout << "s2_at: " << 2*std::abs(s2_at)/n << std::endl;
//   // Rcpp::Rcout << "m2_at: " << 2*std::abs(m2_at)/n << std::endl;
//   // Rcpp::Rcout << "s2_et: " << 2*std::abs(s2_et)/n << std::endl;
//   // Rcpp::Rcout << "m2_et: " << 2*std::abs(m2_et)/n << std::endl;
//   //
//   //
//   // Rcpp::Rcout << "s2_gw: " << s2_gw << std::endl;
//   // Rcpp::Rcout << "m2_gw: " << m2_gw << std::endl;
//   // Rcpp::Rcout << "s2_at: " << s2_at << std::endl;
//   // Rcpp::Rcout << "m2_at: " << m2_at << std::endl;
//   // Rcpp::Rcout << "s2_et: " << s2_et << std::endl;
//   // Rcpp::Rcout << "m2_et: " << m2_et << std::endl;
//
//   return(be);
//
// }
//
//
// // [[Rcpp::export]]
// double be_rau_cpp(Eigen::MatrixXd& x,
//                       const Eigen::VectorXi& spans,
//                       bool detrend,
//                       bool demean,
//                       double taper,
//                       double f1,
//                       double f2,
//                       double frequency_scale) {
//
//
//   if (spans.size() < 1) {
//     Rcpp::stop("spec_pgram: spans must be length 1 or larger.");
//   }
//
//   // detrend or demean
//   x = detrend_and_demean_matrix(x, detrend, demean);
//
//   size_t n_row = x.rows();
//   size_t n_new = next_n_eigen(n_row);
//
//   std::complex<double> scale = 1.0 / n_row; // or n_new
//
//   // taper vector
//   ArrayXd taper_array = spec_taper(n_row, taper).array();
//
//
//   // Do FFTs
//   MatrixXcd pgram = fft_matrix(x.array().colwise() * taper_array,
//                                    n_row);
//
//   unsigned int n = pgram.rows();
//   Eigen::VectorXd freqs = determine_frequency(n) * frequency_scale;
//   Eigen::Vector2i r = get_peaks(freqs, f1, f2);
//
//   // get order of frequencies
//   unsigned int s = 1;
//   unsigned int m = 0;
//   if (r[0] > r[1]) {
//     s = 0;
//     m = 1;
//   }
//
//   std::complex<double> s2_gw, m2_gw, s2_at, m2_at, s2_et, m2_et, d_phase;
//
//   s2_gw = pgram(r[s], 0);
//   m2_gw = pgram(r[m], 0);
//
//   s2_at = pgram(r[s], 1);
//   m2_at = pgram(r[m], 1);
//
//   s2_et = pgram(r[s], 2);
//   m2_et = pgram(r[m], 2);
//
//
//   double be = be_rau_calc_cpp(s2_gw,
//                               s2_et,
//                               s2_at,
//                               m2_gw,
//                               m2_et,
//                               1.0);
//
//
//   // Rcpp::Rcout << "s2_gw: " << s2_gw << std::endl;
//   // Rcpp::Rcout << "m2_gw: " << m2_gw << std::endl;
//   // Rcpp::Rcout << "s2_at: " << s2_at << std::endl;
//   // Rcpp::Rcout << "m2_at: " << m2_at << std::endl;
//   // Rcpp::Rcout << "s2_et: " << s2_et << std::endl;
//   // Rcpp::Rcout << "m2_et: " << m2_et << std::endl;
//
//   return(be);
//
// }




// [[Rcpp::export]]
Rcpp::List be_harmonic_cpp(Eigen::VectorXcd x,
                           bool inverse) {

  // groundwater
  std::complex<double> m2_gw = x[0];
  std::complex<double> s2_gw = x[3];

  // atmospheric pressure
  std::complex<double> m2_at = x[1];
  std::complex<double> s2_at = x[4];

  // earth tides
  std::complex<double> m2_et = x[2];
  std::complex<double> s2_et = x[5];


  double d_phase = std::arg(s2_at) - std::arg(s2_et);
  double term_acworth = 0.0;
  std::complex<double> term_rau(0.0, 0.0);
  double ratio, acworth, rau;

  // Rcpp::Rcout << "phase: " << std::arg(s2_at)-std::arg(s2_gw) << std::endl;
  // Rcpp::Rcout << "dphase: " << d_phase << std::endl;


  //----------------------------------------------------------------------------
  // ratio
  ratio = std::abs(s2_gw) / std::abs(s2_at);


  //----------------------------------------------------------------------------
  // acworth

  if ((std::abs(m2_et) != 0.0) & (std::abs(s2_et) != 0.0)) {

    term_acworth = std::abs(s2_et) * (std::cos(d_phase) * std::abs(m2_gw) / std::abs(m2_et));

    if(inverse) {
      term_acworth = - term_acworth;
    }
  }
  acworth = (std::abs(s2_gw) - term_acworth) / std::abs(s2_at);


  //----------------------------------------------------------------------------
  // rau
  if ((std::abs(m2_et) != 0.0) & (std::abs(s2_et) != 0.0)) {
    term_rau = (m2_gw / m2_et);
  }
  rau = std::abs((s2_gw - term_rau * s2_et) / s2_at);

  //----------------------------------------------------------------------------

  if (inverse) {
    ratio = 1.0 - ratio;
    acworth = 1.0 - acworth;
    rau = 1.0 - rau;
  }

  return(Rcpp::List::create(Named("ratio") = ratio,
                            _["acworth"] = acworth,
                            _["rau"] = rau));



}


//==============================================================================
// [[Rcpp::export]]
Eigen::MatrixXcd be_transfer(Eigen::MatrixXd& x,
                             const Eigen::VectorXi& spans,
                             bool detrend,
                             bool demean,
                             double taper,
                             double frequency,
                             double cycle_size) {

  Eigen::MatrixXcd pgram = spec_pgram(x, spans, detrend, demean, taper);

  unsigned int n = pgram.rows();
  unsigned int frequency_index = std::round(frequency * (n / cycle_size));
  // Rcpp::Rcout << "n: " << n << std::endl;
  // Rcpp::Rcout << "n: " << n << std::endl;
  // Rcpp::Rcout << "cycle_size: " << cycle_size << std::endl;
  // Rcpp::Rcout << "frequency_index: " << frequency_index << std::endl;

  Eigen::MatrixXcd out = solve_cplx_parallel(pgram.row(frequency_index));
  return(out);
}
//==============================================================================

