// #include "hydrorecipes.h"
//
// ==============================================================================
// //' @title
//  //' grf_time
//  //'
//  //' @description
//  //' Parallel convolution of GRF well function and flow rates in the time domain.
//  //' Time series needs to be regularily spaced and so are the flow rates.  Some
//  //' performance gains can be achieved if the number of flow rate does not change
//  //' for each time.
//  //'
//  //' @param radius distance to monitoring interval
//  //' @param storativity aquifer storativity
//  //' @param K aquifer hydraulic conductivity
//  //' @param thickness aquifer thickness
//  //' @param time prediction times
//  //' @param flow_rate well flow rates
//  //' @param flow_time_interval time between flow rate measurements in samples
//  //' @param flow_dimension flow dimension
//  //'
//  //' @return theis solution for multiple pumping scenario
//  //'
//  //'
//  //' @export
//  //'
//  // [[Rcpp::export]]
//  Eigen::VectorXd grf_time(double radius,
//                               double storativity,
//                               double K,
//                               double thickness,
//                               const Eigen::VectorXd &time,
//                               const Eigen::VectorXd &flow_rate,
//                               int flow_time_interval,
//                               double flow_dimension)
//  {
//
//    const double v = (flow_dimension / 2.0) - 1.0;
//
//    size_t n = time.size();
//
//    double coef;
//    Eigen::VectorXd wf(n);
//
//
//    // std::vector<double> u = std::vector<double>(n);
//    // std::vector<double> s = std::vector<double>(n);
//
//    // check that flow time interval values are reasonable
//    if (flow_time_interval < 1)
//    {
//      flow_time_interval = 1;
//      Rcpp::warning("flow_time_interval should be >= 1: set to 1");
//    }
//
//    if (flow_time_interval > n)
//    {
//      flow_time_interval = n;
//      Rcpp::warning("flow_time_interval should be <= n: set to n ");
//    }
//
//
//    for (size_t i = 0; i < n; ++i)
//    {
//      wf[i] = specialfunctions::gamma_inc_single(
//        grf_u(radius, storativity, K, time[i]), v);
//    }
//
//    = grf_coefficient(flow_rate, radius, K, thickness, flow_dimension);
//    // calculate the Well function
//    // u = grf_u_time(radius, storativity, K, time);
//
//    // incomplete gamma function
//    // u = grf_parallel(u, v);
//
//    // calculate the pulse
//    time = impulse_function(time, flow_time_interval);
//
//    // s = well_function_convolve(flow_time_interval, coefs, u);
//
//    return time;
//  }
//
// // //==============================================================================
// // //' @title
// // //' hantush_time_parallel
// // //'
// // //' @description
// // //' Parallel convolution of hantush well function and flow rates in the time domain.
// // //' Time series needs to be regularily spaced.
// // //'
// // //' @param radius distance to monitoring interval
// // //' @param storativity aquifer storativity
// // //' @param transmissivity aquifer transmissivity
// // //' @param leakage hantush leakage
// // //' @param time prediction times
// // //' @param flow_rate well flow rates
// // //' @param flow_time_interval time between flow rate measurements in samples
// // //' @param n_terms number of terms to use in Hantush solution.  More is more precise but slower.
// // //'
// // //' @return hantush jacob solution for multiple pumping scenario
// // //'
// // //'
// // //' @export
// // //'
// // // [[Rcpp::export]]
// // Rcpp::NumericVector hantush_time_parallel(double radius,
// //                                           double storativity,
// //                                           double transmissivity,
// //                                           double leakage,
// //                                           const Rcpp::NumericVector& time,
// //                                           const Rcpp::NumericVector& flow_rate,
// //                                           int flow_time_interval,
// //                                           int n_terms) {
//
// //   int n = time.size();
// //   double b = hantush_epsilon(radius, leakage);
//
// //   const Rcpp::NumericVector coefs = well_function_coefficient(flow_rate,  transmissivity);
//
// //   Rcpp::NumericVector u = Rcpp::NumericVector(n);
// //   Rcpp::NumericVector s = Rcpp::NumericVector(n);
//
// //   // calculate the Well function
// //   u = grf_u_time(radius, storativity, transmissivity, time);
// //   u = hantush_well_parallel(u, b, n_terms);
//
// //   // calculate the pulse
// //   u = impulse_function(u, flow_time_interval);
//
// //   s = well_function_convolve(flow_time_interval, coefs, u);
//
// //   return s;
// // }
//
// // //==============================================================================
// // //' @title
// // //' grf
// // //'
// // //' @description
// // //' Parallel convolution of GRF well function and flow rates in the time domain.
// // //' Time series needs to be regularily spaced and so are the flow rates.  Some
// // //' performance gains can be achieved if the number of flow rate does not change
// // //' for each time.
// // //'
// // //' @param radius distance to monitoring interval
// // //' @param storativity aquifer storativity
// // //' @param transmissivity aquifer transmissivity
// // //' @param leakage hantush leakage
// // //' @param time prediction times
// // //' @param flow_rate well flow rates
// // //' @param flow_rate_times time of flow rate changes
// // //' @param n_terms number of terms to use in Hantush solution.  More is more precise but slower.
// // //'
// // //' @return theis solution for multiple pumping scenario
// // //'
// // //'
// // //' @export
// // //'
// // // [[Rcpp::export]]
// // Rcpp::NumericVector hantush(double radius,
// //                             double storativity,
// //                             double transmissivity,
// //                             double leakage,
// //                             const Rcpp::NumericVector& time,
// //                             const Rcpp::NumericVector& flow_rate,
// //                             const Rcpp::NumericVector& flow_rate_times,
// //                             int n_terms) {
//
// //   int n = time.size();
// //   double b = hantush_epsilon(radius, leakage);
// //   int n_flow_rate = flow_rate.size();
//
// //   const Rcpp::NumericVector coefs = well_function_coefficient(flow_rate,  transmissivity);
//
// //   Rcpp::NumericVector u = Rcpp::NumericVector(n);
// //   Rcpp::NumericVector s = Rcpp::NumericVector(n);
// //   Rcpp::NumericVector t_elapsed = Rcpp::NumericVector(n);
// //   Rcpp::LogicalVector idx = Rcpp::LogicalVector(n);
//
// //   for (int j = 0; j < n_flow_rate; j++) {
//
// //     // get the elapsed time
// //     t_elapsed = time - flow_rate_times[j];
// //     idx = t_elapsed > 0;
//
// //     // calculate the Well function
// //     u[idx] = grf_u_time(radius, storativity, transmissivity, t_elapsed[idx]);
//
// //     // incomplete gamma function
// //     u[idx] = hantush_well_parallel(u[idx], b, n_terms);
// //     u[!idx] = 0;
//
// //     //impulse
// //     if (j > 0) {
// //       s = s + (coefs[j] * u) - (coefs[j-1] * u);
// //     } else {
// //       s = (coefs[j] * u);
// //     }
//
// //   }
//
// //   return s;
// // }
//
//
//
// // //==============================================================================
// // //' @title
// // //' grf
// // //'
// // //' @description
// // //' Parallel convolution of GRF well function and flow rates in the time domain.
// // //' Time series needs to be regularily spaced and so are the flow rates.  Some
// // //' performance gains can be achieved if the number of flow rate does not change
// // //' for each time.
// // //'
// // //' @param radius distance to monitoring interval
// // //' @param storativity aquifer storativity
// // //' @param K aquifer hydraulic conductivity
// // //' @param thickness aquifer thickness
// // //' @param time prediction times
// // //' @param flow_rate well flow rates
// // //' @param flow_rate_times times where flow rates change
// // //' @param flow_dimension flow dimension
// // //'
// // //' @return theis solution for multiple pumping scenario
// // //'
// // //'
// // //' @export
// // //'
// // // [[Rcpp::export]]
// // Rcpp::NumericVector grf(double radius,
// //                         double storativity,
// //                         double K,
// //                         double thickness,
// //                         const Rcpp::NumericVector& time,
// //                         const Rcpp::NumericVector& flow_rate,
// //                         const Rcpp::NumericVector& flow_rate_times,
// //                         double flow_dimension) {
//
// //   double v = (flow_dimension / 2.0) - 1.0;
//
// //   int n = time.size();
// //   int n_flow_rate = flow_rate.size();
//
// //   const Rcpp::NumericVector coefs = grf_coefficient(flow_rate, radius, K, thickness, flow_dimension);
//
// //   Rcpp::NumericVector u = Rcpp::NumericVector(n);
// //   Rcpp::NumericVector s = Rcpp::NumericVector(n);
// //   Rcpp::NumericVector t_elapsed = Rcpp::NumericVector(n);
// //   Rcpp::LogicalVector idx = Rcpp::LogicalVector(n);
//
// //   for (int j = 0; j < n_flow_rate; j++) {
//
// //     // get the elapsed time
// //     t_elapsed = time - flow_rate_times[j];
// //     idx = t_elapsed > 0;
//
// //     // calculate the Well function
// //     u[idx] = grf_u_time(radius, storativity, K, t_elapsed[idx]);
//
// //     // incomplete gamma function
// //     u[idx] = grf_parallel(u[idx], v);
// //     u[!idx] = 0;
//
// //     //impulse
// //     if (j > 0) {
// //       s = s + (coefs[j] * u) - (coefs[j-1] * u);
// //     } else {
// //       s = (coefs[j] * u);
// //     }
//
// //   }
//
// //   return s;
// // }
//
//
// // // FFTW is about 3 times faster-------------------------------------------------
//
// // // //==============================================================================
// // // //' @title
// // // //' hantush_freq
// // // //'
// // // //' @description
// // // //' Parallel convolution of hantush well function and flow rates in the frequency domain.
// // // //' Time series needs to be regularily spaced.
// // // //'
// // // //' @param radius distance to monitoring interval
// // // //' @param storativity aquifer storativity
// // // //' @param transmissivity aquifer transmissivity
// // // //' @param leakage hantush leakage
// // // //' @param time prediction times
// // // //' @param flow_rate well flow rates
// // // //' @param n_terms number of terms to use in Hantush solution.  More is more precise but slower.
// // // //'
// // // //' @return hantush jacob solution for multiple pumping scenario
// // // //'
// // // //'
// // // //' @export
// // // //'
// // // // [[Rcpp::export]]
// // // arma::colvec hantush_freq(double radius,
// // //                           double storativity,
// // //                           double transmissivity,
// // //                           double leakage,
// // //                           const Rcpp::NumericVector& time,
// // //                           const Rcpp::NumericVector& flow_rate,
// // //                           int n_terms) {
// // //
// // //   int n = time.size();
// // //   int n_q = flow_rate.size();
// // //   double b = hantush_epsilon(radius, leakage);
// // //
// // //   Rcpp::NumericVector coefs = well_function_coefficient(flow_rate,  transmissivity);
// // //
// // //   Rcpp::NumericVector u = Rcpp::NumericVector(n);
// // //   //Rcpp::NumericVector s = Rcpp::NumericVector(n);
// // //   arma::colvec s(n);
// // //
// // //   // calculate the Well function
// // //   u = theis_u_time(radius, storativity, transmissivity, time);
// // //   u = hantush_well_parallel(u, b, n_terms);
// // //
// // //   // calculate the pulse
// // //   u = impulse_function(u, 1);
// // //
// // //   arma::colvec y(u.begin(), u.size(), false);
// // //   arma::colvec z(coefs.begin(), coefs.size(), false);
// // //
// // //   s = arma::real(arma::ifft(arma::fft(z) % arma::fft(y)));
// // //
// // //   return s;
// // // }
// // //
// // //
// // // //==============================================================================
// // // //' @title
// // // //' grf_freq
// // // //'
// // // //' @description
// // // //' Parallel convolution of Theis well function and flow rates in the time domain.
// // // //' Time series needs to be regularily spaced.
// // // //'
// // // //' @param radius distance to monitoring interval
// // // //' @param storativity aquifer storativity
// // // //' @param transmissivity aquifer transmissivity
// // // //' @param time prediction times
// // // //' @param flow_rate well flow rates
// // // //' @param flow_dimension time between flow rate measurements in samples
// // // //'
// // // //' @return theis solution for multiple pumping scenario
// // // //'
// // // //'
// // // //' @export
// // // //'
// // // // [[Rcpp::export]]
// // // arma::colvec grf_freq(double radius,
// // //                       double storativity,
// // //                       double transmissivity,
// // //                       const Rcpp::NumericVector& time,
// // //                       const Rcpp::NumericVector& flow_rate,
// // //                       double flow_dimension) {
// // //
// // //   double v = 1 - flow_dimension/2;
// // //
// // //   int n = time.size();
// // //   int n_q = flow_rate.size();
// // //
// // //   Rcpp::NumericVector coefs = well_function_coefficient(flow_rate,  transmissivity);
// // //
// // //   Rcpp::NumericVector u = Rcpp::NumericVector(n);
// // //   arma::colvec s(n);
// // //
// // //   // calculate the Well function
// // //   u = theis_u_time(radius, storativity, transmissivity, time);
// // //   u = grf_parallel(v, u);
// // //
// // //   // calculate the pulse
// // //   u = impulse_function(u, 1);
// // //
// // //   arma::colvec y(u.begin(), u.size(), false);
// // //   arma::colvec z(coefs.begin(), coefs.size(), false);
// // //
// // //   s = arma::real(arma::ifft(arma::fft(z) % arma::fft(y)));
// // //
// // //   return s;
// // // }
