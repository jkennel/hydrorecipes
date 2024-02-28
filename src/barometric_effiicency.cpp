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

  int n = dep.n_elem-lag_space;
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
//' be_acworth_cpp need to work on this
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
double be_acworth_cpp(const double s2_gw,
                      const double s2_et,
                      const double s2_at,
                      const double m2_gw,
                      const double m2_et,
                      const double d_phase,
                      const bool inverse) {

  const double term = s2_et * (cos(d_phase) * m2_gw / m2_et);

  if (inverse) {
    return ((s2_gw + term) / s2_at);
  }

  return ((s2_gw - term) / s2_at);
}
