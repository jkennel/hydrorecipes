#include "hydrorecipes.h"

// Ogata, A., Banks, R.B., 1961. A solution of the differential equation of
// longitudinal dispersion in porous media. U. S. Geol. Surv. Prof. Pap. 411-A.
//
// 1-D
// infinite source
// uniform flow
// constant parameters
// no decay
// no retardation
//==============================================================================
//' @title
//' Ogata-Banks solution for 1-D flow.
//'
//' @description
//' Ogata, A., Banks, R.B., 1961. A solution of the differential equation of
//' longitudinal dispersion in porous media. U. S. Geol. Surv. Prof. Pap. 411-A.
//' 1-D, infinite source, uniform flow, constant parameters, no decay, no retardation
//'
//' @param D diffusion coefficient
//' @param v double velocity
//' @param C0 double concentration
//' @param x double x position
//' @param t double time
//'
//' @return ogata banks solution
//'
//' @export
//'
//' @noRd
//'
// [[Rcpp::export]]
double ogata_banks_ind(const double D,
                       const double v,
                       const double C0,
                       const double x,
                       double t) {

  const double term_1 = erfc((x - v * t) / (2 * sqrt(D * t)));
  const double term_2 = erfc((x + v * t) / (2 * sqrt(D * t)));

  if (term_2 == 0.0) {
    return(0.5 * C0 * term_1);
  }

  const double term_3 = exp(v * x / D);

  return 0.5 * C0 * (term_1 + term_3 * term_2);

}


// [[Rcpp::export]]
std::vector<double> ogata_banks_vec(const double D,
                                    const double v,
                                    const double C0,
                                    const double x,
                                    std::vector<double> t) {

  for (auto &out : t)
    out = ogata_banks_ind(D, v, C0, x, out);

  return(t);

}


// [[Rcpp::export]]
double ogata_banks_decay_ind(
    const double c0,
    const double v,
    const double D,
    const double R,
    const double k,
    const double x,
    const double t
) {

  if(t <= 0.0) {
    return(0.0);
  }

  const double B =    sqrt(pow((v / (2.0 * D)), 2) + (k * R / D));
  const double term = sqrt(pow((v / R), 2) + (4.0 * k * D / R));

  const double exp_x = std::exp(B * x);
  const double exp_rx = 1.0 / exp_x;

  const double term_t = term * t;
  const double co = 1.0 / (2.0 * sqrt(D * t / R));


  double erfc_1 = 0.0;
  double erfc_2 = 0.0;
  double output = 0.0;

  // check for overflow
  if ((exp_rx != 0.0) & ((x - term_t) < 30)) {
    erfc_1 = (exp_rx * std::erfc((x - term_t) * co));
  }

  if ((exp_x != 0.0) & ((x + term_t) < 30)) {
    erfc_2 = (exp_x  * std::erfc((x + term_t) * co));
  }

  if(std::isnan(erfc_1)) {
    erfc_1 = 0.0;
  }

  if(std::isnan(erfc_2)) {
    erfc_2 = 0.0;
  }

  if((erfc_1 != 0.0) || (erfc_2 != 0.0)) {
    output = 0.5 * c0 * exp(v * x / (2.0 * D)) * (erfc_1 + erfc_2);
  }

  // Rcpp::Rcout << "The value B " << B << std::endl;
  // Rcpp::Rcout << "The value co " << co << std::endl;
  // Rcpp::Rcout << "The value exp_rx " << exp_rx << std::endl;
  // Rcpp::Rcout << "The value exp_x " << exp_x << std::endl;
  // Rcpp::Rcout << "The value term_t " << term_t << std::endl;
  // Rcpp::Rcout << "The value x-term_t " << x-term_t << std::endl;

  // Rcpp::Rcout << "The value exp(v * x / (2.0 * D))" << exp(v * x / (2.0 * D)) << std::endl;
  // Rcpp::Rcout << "The value (exp_rx * std::erfc((x - term_t) * co))" << (exp_rx * std::erfc((x - term_t) * co)) << std::endl;
  // Rcpp::Rcout << "The value (exp_x  * std::erfc((x + term_t) * co))" << (exp_x  * std::erfc((x + term_t) * co)) << std::endl;
  // Rcpp::Rcout << "The value erfc_1" << erfc_1 << std::endl;
  // Rcpp::Rcout << "The value erfc_2" << erfc_2 << std::endl;


  return output;

}

// [[Rcpp::export]]
Rcpp::NumericVector ogata_banks_decay_vec(
    const double c0,
    const double v,
    const double D,
    const double R,
    const double k,
    Rcpp::NumericVector x,
    Rcpp::NumericVector t
) {

  unsigned int n = x.size();
  Rcpp::NumericVector out(n);

  for (unsigned int i = 0; i < n; ++i) {
    out[i] = ogata_banks_decay_ind(c0, v, D, R, k, x[i], t[i]);
  }

  return(out);
}


// // Ogata, A., Banks, R.B., 1961. A solution of the differential equation of
// // longitudinal dispersion in porous media. U. S. Geol. Surv. Prof. Pap. 411-A.
// //
// // 1-D
// // infinite source
// // uniform flow
// // constant parameters
// // no decay
// // no retardation
// //==============================================================================
// //' @title
// //' Ogata-Banks solution for 1-D flow (vectorized).
// //'
// //' @description
// //' Ogata, A., Banks, R.B., 1961. A solution of the differential equation of
// //' longitudinal dispersion in porous media. U. S. Geol. Surv. Prof. Pap. 411-A.
// //' 1-D, infinite source, uniform flow, constant parameters, decay, retardation
// //'
// //' To have values match the excel sheet
// //' https://www.civil.uwaterloo.ca/jrcraig/pdf/OgataBanks.xlsm the decay
// //' coefficient needs to be scaled by the retardation coefficient.
// //'
// //' @param diffusion double diffusion coefficient
// //' @param retardation double retardation coefficient
// //' @param decay double decay coefficient
// //' @param velocity double velocity
// //' @param concentration_initial double concentration
// //' @param distance vector x position
// //' @param time vector time
// //'
// //' @return ogata banks solution each row is an x value and each column is a time
// //'
// //' @export
// //'
// // [[Rcpp::export]]
// Rcpp::List ogata_banks_list(
//     Eigen::ArrayXd time,
//     Eigen::ArrayXd distance,
//     double concentration_initial,
//     double velocity,
//     double diffusion,
//     double retardation,
//     double decay
// ) {
//
//
//   const double B = sqrt(pow((velocity / (2.0 * diffusion)), 2) +
//                         (decay * retardation / diffusion));
//   const double term = sqrt(pow((velocity / retardation), 2) + (4.0 * decay * diffusion / retardation));
//
//   const Eigen::ArrayXd exp_x = (B * distance).exp();
//   const Eigen::ArrayXd exp_rx = 1.0 / exp_x;
//   const Eigen::ArrayXd term_t = term * time;
//   const Eigen::ArrayXd co = 1.0 / (2.0 * (diffusion * time / retardation).sqrt());
//
//   // this should be changed to just make one column
//   // for (unsigned int i = 0; i < n_t; ++i) {
//   const Eigen::ArrayXd output = 0.5 * concentration_initial * (velocity * distance / (2.0 * diffusion)).exp() *
//     ((exp_rx * ((distance - term_t) * co).erfc()) +
//     (exp_x * ((distance + term_t) * co).erfc()));
//
//   // }
//   return Rcpp::List::create(Rcpp::Named("ogata_banks") = output);
//
//   // return(output);
//
// }

//
// Eigen::MatrixXd ogata_banks(double D, double R, double decay,
//                              double v, double C0,
//                              Eigen::ArrayXd x,
//                              Eigen::ArrayXd t) {
//
//   const unsigned int n_x = x.size();
//   const unsigned int n_t = t.size();
//
//   const double B = sqrt(pow((v / (2.0 * D)), 2) + (decay * R / D));
//   const double term = sqrt(pow((v / R), 2) + (4.0 * decay * D / R));
//
//   Eigen::ArrayXXd output(n_x, n_t);
//   Eigen::ArrayXXd output1(n_x, n_t);
//   Eigen::ArrayXXd output2(n_x, n_t);
//
//   const Eigen::ArrayXd exp_x = (B * x).exp();
//   const Eigen::ArrayXd exp_rx = 1.0 / exp_x;
//
//   const Eigen::RowVectorXd term_t = term * t;
//   const Eigen::RowVectorXd co = 1.0 / (2.0 * (D * t / R).array().sqrt());
//
//   output.colwise() = x;
//   output1.colwise() = exp_x;
//   output2.colwise() = exp_rx;
//
//   output = 0.5 * C0 * (v * output / (2.0 * D)).exp() *
//
//     ((output2 *
//     ((output.rowwise() - term_t.array()).rowwise() * co.array()).erfc()) +
//
//     (output1 *
//     ((output.rowwise() + term_t.array()).rowwise() * co.array()).erfc()));
//
//   return(output.matrix());
//
// }
//
// // [[Rcpp::export]]
// Rcpp::NumericVector pp(Rcpp::NumericVector x) {
//
//   return(2.0 * Rcpp::pnorm5(x * sqrt(2.0), 0, 1, false, false));
// }
//
//
//
//
//
//  // [[Rcpp::export]]
//  arma::mat ogata_banks3(double D, double R, double decay,
//                        double v, double C0,
//                        arma::vec x,
//                        arma::rowvec t) {
//
//    int n_x = x.n_elem;
//    int n_t = t.n_elem;
//
//    double B = sqrt(pow((v / (2.0 * D)), 2) + (decay * R / D));
//    double term = sqrt(pow((v / R), 2) + (4.0 * decay * D / R));
//
//    arma::mat output(n_x, n_t);
//
//    for (int i = 0; i < n_t; i++) {
//      output.col(i) = x;
//    }
//
//    output = 0.5 * C0 * exp(v * output / (2.0 * D)) %
//      ((exp(-B * output) %
//      arma::erfc((output.each_row() - term * t).each_row() / (2.0 * sqrt(D * t / R)))) +
//      (exp(B * output) %
//      arma::erfc((output.each_row() + term * t).each_row() / (2.0 * sqrt(D * t / R)))));
//
//    return(output);
//
//  }
//
// // [[Rcpp::export]]
// Eigen::VectorXd eer(Eigen::VectorXd x) {
//   return(x.array().erfc());
// }
//
// // [[Rcpp::export]]
// std::vector<double> er(std::vector<double> x) {
//   return(specialfunctions::erfc_vec(x));
// }
//
// // [[Rcpp::export]]
// arma::vec aer(arma::vec x) {
//
//   return(arma::erfc(x));
//
// }
//
//
// // [[Rcpp::export]]
// arma::mat ogata_banks2(double D, double R, double decay,
//                        double v, double C0,
//                        arma::vec x,
//                        arma::rowvec t) {
//
//   int n_x = x.n_elem;
//   int n_t = t.n_elem;
//
//   double B = sqrt(pow((v / (2.0 * D)), 2) + (decay * R / D));
//   double term = sqrt(pow((v / R), 2) + (4.0 * decay * D / R));
//
//   arma::mat output(n_x, n_t);
//   arma::mat output1(n_x, n_t);
//   arma::mat output2(n_x, n_t);
//
//   arma::vec exp_x = exp(B * x);
//   arma::vec exp_rx = 1.0 / exp_x;
//
//   arma::rowvec term_t = term * t;
//   arma::rowvec co = (2.0 * sqrt(D * t / R));
//   for (int i = 0; i < n_t; i++) {
//     output.col(i) = x;
//     output1.col(i) = exp_x;
//     output2.col(i) = exp_rx;
//   }
//
//   output = 0.5 * C0 * exp(v * output / (2.0 * D)) %
//     ((output2 %
//     arma::erfc((output.each_row() - term_t).each_row() / co)) +
//     (output1 %
//     arma::erfc((output.each_row() + term_t).each_row() / co)));
//
//   return(output);
//
// }
//
//
//
//
// // [[Rcpp::export]]
// Eigen::MatrixXd vop(Eigen::VectorXd x, Eigen::RowVectorXd y) {
//   return (x * y);
// }

/*** R
x <- abs(rnorm(2000000))
bench::mark(
  pracma::erfc(x),
  hydrorecipes:::eer(x),
  hydrorecipes:::er(x),
  (hydrorecipes:::aer(x)),
  hydrorecipes:::pp(x),
  check = FALSE
)

a <- rnorm(10000000)

bench::mark(
  hydrorecipes:::pp(a),
  hydrorecipes:::eer(a),
  check = FALSE

)

D <- 0.1
R <- 1
decay <- 0
v <- 0.1
C0 <- 1
x <- 1:50000
t <- 1:50

bench::mark(
  # hydrorecipes:::ogata_banks(D = D, R = R, decay = decay,
  #                        v = v, C0 = C0, x = x,
  #                        t = t),
  # hydrorecipes:::ogata_banks2(D = D, R = R, decay = decay,
  #                         v = v, C0 = C0, x = x,
  #                         t = t),
  # hydrorecipes:::ogata_banks3(D = D, R = R, decay = decay,
  #                         v = v, C0 = C0, x = x,
  #                         t = t),
  hydrorecipes:::ogata_banks_list(D = D, R = R, decay = decay,
                              v = v, C0 = C0, x = x,
                              t = t),
  check = FALSE
  # relative = TRUE
)


*/

