#include "hydrorecipes.h"

// Bouwer, H., 1989. The Bouwer and Rice Slug Test—An Update a. Groundwater, 27(3), pp.304-309.

// http://hydrotools.sourceforge.net/
// [[Rcpp::export]]
double a_cpp(double x) {

  if(x < 2.554422663){
    x = 1.638445671 + 0.166908063 * x + 0.000740459 *
      std::exp(6.17105281 * x - 1.054747686 * x * x);
  } else {
    x = 11.00393028 - 170.7752217 * std::exp(-1.509639982 * x);
  }

  return(x);

}

// http://hydrotools.sourceforge.net/
// [[Rcpp::export]]
double b_cpp(double x) {

  if(x < 2.596774459){
    x = 0.174811819 + 0.060059188 * x + 0.007965502 *
      std::exp(2.053376868 * x - 0.007790328 * x * x);
  } else {
    x = 4.133124586 - 93.06136936 * std::exp(-1.435370997 * x);
  }


  return(x);

}

// http://hydrotools.sourceforge.net/
// [[Rcpp::export]]
double c_cpp(double x) {

  if(x < 2.200426117){
    x = 0.074711376 + 1.083958569 * x + 0.00557352 *
      std::exp(2.929493814 * x - 0.001028433 * x * x);

  } else {
    x = 15.66887372 - 178.4329289 * std::exp(-1.322779744 * x);
  }

  return(x);
}


//==============================================================================
//' Calculate equations 4 and 5 from bouwer, 1989
//'
//' @param rw radius of well
//' @param Le Effecive screen length
//' @param Lw height of water from bottom of well
//' @param H height from bottom of aquifer
//'
//' @return ln(Re/rw)
//'
//' @noRd
// [[Rcpp::export]]
double bouwer_rice_abc(double rw, double Le, double Lw, double H) {

  double lr, llr, a, b, c, ret, t1, t2;

  lr  = Le / rw;
  llr = std::log10(lr);

  // calculate empirical values a, b, c
  a = a_cpp(llr);
  b = b_cpp(llr);
  c = c_cpp(llr);

  Rcpp::Rcout << "a: " << a << std::endl;
  Rcpp::Rcout << "b: " << b << std::endl;
  Rcpp::Rcout << "c: " << c << std::endl;

  t1 = (1.1 / std::log(Lw / rw));
  t2 = std::log((H - Lw) / rw);

  // calculate ln(Re/rw)
  if (Lw < H){
    ret = std::pow(t1 + (a + b * t2)/lr, -1);
  } else {
    ret = std::pow(t1 + c/lr, -1);
  }

  return(ret);
}


//==============================================================================
//' Calculate transmissivity with Bouwer-Rice solution
//'
//' @param time the elapsed time
//' @param drawdown the drawdown
//' @param radius_screen radius of the screen
//' @param radius_casing radius of the casing where the water level is
//' @param Le Effecive screen length
//' @param Lw height of water from bottom of well
//' @param H height from bottom of aquifer
//'
//' @return transmissivity from bouwer_rice
//'
//' @export
//'
//' @noRd
// [[Rcpp::export]]
double bouwer_rice(arma::vec time,
                   arma::vec drawdown,
                   double radius_screen,
                   double radius_casing,
                   double Le,
                   double Lw,
                   double H) {

  arma::mat x(time.n_elem, 2);
  x.fill(1.0);
  x.col(0) = time;

  double wc = bouwer_rice_abc(radius_screen, Le, Lw, H);
  double drawdown0  = drawdown(0);

  arma::vec slope = arma::solve(x, arma::log10(drawdown / drawdown0));

  // Calculate the hydraulic conductivity
  return((radius_casing * radius_casing * wc) * (-slope[0] * 2.302585) /
         (2.0 * Le) );

}


/***R
# library(data.table)
# data(bouwer)
# bw <- as.data.table(bouwer)
#
# rc <- 4/2/12        # radius of 2 inches
# rw <- 8.25/2/12     # radius of screen
# Le <- 10            # screen length
# y0 <- 1.44          # initial drawdown
# Lw <- 17.92         # height of water above screen bottom
# wl <- 6.08          # static wl
# H  <- 18.92         # height of water level above base
# t  <- as.numeric(bouwer$datetime - bouwer$datetime[1])    # elapsed time
# y  <- wl - bouwer$val                                     # change in wl
#
#
# a <- bouwer_rice(t, y, rw, rc, Le, Lw, H)
#
# a * 86400

*/
