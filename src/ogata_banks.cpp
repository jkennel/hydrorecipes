#include "frecipes.h"

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
 // [[Rcpp::export]]
 double ogata_banks_ind(double D, double v, double C0, double x,
                        double t) {

   return 0.5 * C0 * (erfc((x - v * t) / (2 * sqrt(D * t))) +
                      exp(v * x / D) * erfc((x + v * t) / (2 * sqrt(D * t))));

 }


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
 //' Ogata-Banks solution for 1-D flow (vectorized).
 //'
 //' @description
 //' Ogata, A., Banks, R.B., 1961. A solution of the differential equation of
 //' longitudinal dispersion in porous media. U. S. Geol. Surv. Prof. Pap. 411-A.
 //' 1-D, infinite source, uniform flow, constant parameters, decay, retardation
 //'
 //' To have values match the excel sheet
 //' https://www.civil.uwaterloo.ca/jrcraig/pdf/OgataBanks.xlsm the decay
 //' coefficient needs to be scaled by the retardation coefficient.
 //'
 //' @param diffusion double diffusion coefficient
 //' @param retardation double retardation coefficient
 //' @param decay double decay coefficient
 //' @param velocity double velocity
 //' @param concentration_initial double concentration
 //' @param distance vector x position
 //' @param time vector time
 //'
 //' @return ogata banks solution each row is an x value and each column is a time
 //'
 //' @export
 //'
 // [[Rcpp::export]]
 Rcpp::List ogata_banks_list(
     Eigen::ArrayXd time,
     Eigen::ArrayXd distance,
     double concentration_initial,
     double velocity,
     double diffusion,
     double retardation,
     double decay
 ) {

   const unsigned int n_x = distance.size();
   const unsigned int n_t = time.size();

   const double B = sqrt(pow((velocity / (2.0 * diffusion)), 2) +
                         (decay * retardation / diffusion));
   const double term = sqrt(pow((velocity / retardation), 2) + (4.0 * decay * diffusion / retardation));

   Rcpp::List output(n_t);

   const Eigen::ArrayXd exp_x = (B * distance).exp();
   const Eigen::ArrayXd exp_rx = 1.0 / exp_x;
   const Eigen::ArrayXd term_t = term * time;
   const Eigen::ArrayXd co = 1.0 / (2.0 * (diffusion * time / retardation).sqrt());


   for (unsigned int i = 0; i < n_t; ++i) {
     output(i) = 0.5 * concentration_initial * (velocity * distance / (2.0 * diffusion)).exp() *
       ((exp_rx * ((distance - term_t(i)) * co(i)).erfc()) +
       (exp_x * ((distance + term_t(i)) * co(i)).erfc()));

   }

   return(output);

 }

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
  frecipes:::eer(x),
  frecipes:::er(x),
  (frecipes:::aer(x)),
  frecipes:::pp(x),
  check = FALSE
)

a <- rnorm(10000000)

bench::mark(
  frecipes:::pp(a),
  frecipes:::eer(a),
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
  frecipes:::ogata_banks(D = D, R = R, decay = decay,
                         v = v, C0 = C0, x = x,
                         t = t),
  # frecipes:::ogata_banks2(D = D, R = R, decay = decay,
  #                         v = v, C0 = C0, x = x,
  #                         t = t),
  # frecipes:::ogata_banks3(D = D, R = R, decay = decay,
  #                         v = v, C0 = C0, x = x,
  #                         t = t),
  frecipes:::ogata_banks_list(D = D, R = R, decay = decay,
                              v = v, C0 = C0, x = x,
                              t = t),
  check = FALSE
  # relative = TRUE
)


*/

