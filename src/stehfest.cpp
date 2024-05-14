// #include "hydrorecipes.h"
//
// using namespace boost::math;
//
// //' stehfest_v
// //'
// //' @description Inverse Laplace transform using Stehfest algorithm
// //'
// //' @param n coefficient for Stehfest algorithm (should be < 20)
// //'
// //' @return inverse laplace transform of impulse function
// //'
// //' @noRd
// //'
// // [[Rcpp::export]]
// Eigen::RowVectorXd stehfest_v(
//      const int n
// ) {
//
//   assert(n < 20);
//   assert(n % 2 == 0);
//
//   int n_div_2 = n / 2;
//   int s, e;
//   double z;
//
//   Eigen::RowVectorXd fact(n+1);
//   Eigen::RowVectorXd v(n);
//
//   // calculate factorials
//   for (unsigned int i = 0; i <= n; ++i) {
//     fact[i] = factorial<double>(i);
//   }
//
//   for (int i = 0; i < n; ++i) {
//     if (i > n_div_2 - 1) {
//       s = n_div_2 - 1;
//     } else {
//       s = i;
//     }
//     e = (int)(i / 2);
//
//     z = 0.0;
//     for (int k = s; k >= e; --k) {
//       z += (pow((double)k + 1.0, n_div_2) * fact[2 * k + 2]) /
//         (fact[n_div_2 - k - 1] * fact[k + 1] * fact[k] *
//           fact[i - k] * fact[2 * k - i + 1]);
//
//     }
//
//     v[i] = pow(-1.0, n_div_2 + i + 1.0) * z;
//
//   }
//
//   return(v);
//
// }
//
//
// // [[Rcpp::export]]
// Eigen::MatrixXd stehfest_p(
//     Eigen::VectorXd time,
//     int n_stehfest
// ) {
//
//   Eigen::RowVectorXd v = Eigen::RowVectorXd::LinSpaced(n_stehfest, 1.0, (double)n_stehfest);
//
//   return(time * v);
//
// }
//
//
// // [[Rcpp::export]]
// double barker_herbert_laplace(double p,
//                               double radius,
//                               double radius_patch,
//                               double t_1,
//                               double t_2,
//                               double s_1,
//                               double s_2) {
//
//   double N = sqrt(s_1 * p / t_1);
//   double A = sqrt(s_2 * p / t_2);
//
//   double ct = (t_2 / t_1) * (A / N);
//
//   double bi_n0 = cyl_bessel_i(0.0, N * radius_patch);
//   double bk_a0 = cyl_bessel_k(0.0, A * radius_patch);
//   double bk_n0 = cyl_bessel_k(0.0, N * radius_patch);
//   double bi_n1 = cyl_bessel_i(1.0, N * radius_patch);
//   double bk_a1 = cyl_bessel_k(1.0, A * radius_patch);
//   double bk_n1 = cyl_bessel_k(1.0, N * radius_patch);
//
//   double denom = (ct * bi_n0 * bk_a1 + bi_n1 * bk_a0) * p;
//
//   double term_1 = (bk_n1 * bk_a0  - bk_a1 * bk_n0 * ct);
//   double term_2 = (bk_n0 * bi_n1 + bk_n1 * bi_n0);
//
//   double drawdown;
//
//   // if radius is inside the patch
//   if (radius <= radius_patch) {
//
//     drawdown = cyl_bessel_k(0.0, N * radius) / p  +
//      (term_1 * cyl_bessel_i(0.0, N * radius)) / denom;
//
//     return(drawdown);
//
//   }
//
//   // if radius is outside the patch
//   drawdown = term_2 * cyl_bessel_k(0.0, A * radius) / denom;
//
//   return(drawdown);
//
// }
//
// // [[Rcpp::export]]
// double hantush_laplace(double p,
//                        double radius,
//                        double radius_patch,
//                        double t_1,
//                        double s_1,
//                        double c,
//                        double q) {
//
//   double w = (s_1 * c * p + 1.0) / (c * t_1);
//   return -q / (2.0 * M_PI * t_1 * p) * cyl_bessel_k(0.0, radius * sqrt(w));
// }
//
// // [[Rcpp::export]]
// double theis_large_diameter_laplace(double p,
//                                     double radius,
//                                     double radius_well,
//                                     double t_1,
//                                     double s_1
// ) {
//
//   double q = sqrt(s_1 * p / t_1);
//   double radius_casing = radius_well;
//
//   double bk_0 = cyl_bessel_k(0.0, q * radius);
//   double bk_rw_0 = cyl_bessel_k(0.0, q * radius_well);
//   double bk_rw_1 = cyl_bessel_k(1.0, q * radius_well);
//
//   double denom = M_PI * p * ((radius_casing * radius_casing * p * bk_rw_0) +
//     (2.0 * radius_well * t_1 * q * bk_rw_1));
//
//
//   return(bk_0 / denom);
//
// }
//
//
// // [[Rcpp::export]]
// Rcpp::List stehfest_barker_herbert(
//      Eigen::VectorXd time,
//      double flow_rate,
//      double radius,
//      double radius_patch,
//      double t_1,
//      double t_2,
//      double s_1,
//      double s_2,
//      unsigned int n_stehfest
// ) {
//
//   // time to tau
//   time = std::log(2.0) / time.array();
//
//   Eigen::VectorXd v = stehfest_v(n_stehfest);
//   Eigen::MatrixXd p = stehfest_p(time, n_stehfest);
//
//   std::vector<double> p_vec(p.data(), p.data() + p.size());
//
//   for (auto &out : p_vec)
//     out = barker_herbert_laplace(out, radius, radius_patch, t_1, t_2, s_1, s_2);
//
//   p = Eigen::Map<Eigen::MatrixXd>(p_vec.data(), p.rows(), p.cols());
//
//   p.array() *= flow_rate / (2.0 * M_PI * t_1);
//
//   Eigen::VectorXd ret = (p * v).array() * time.array();
//
//   return Rcpp::List::create(Rcpp::Named("barker_herbert") = ret);
//
// }
//
//
//
//
//
//
//
//
//
// // [[Rcpp::export]]
// Rcpp::List stehfest_theis_large_diameter(
//      Eigen::VectorXd time,
//      double flow_rate,
//      double radius,
//      double radius_well,
//      double t_1,
//      double s_1,
//      unsigned int n_stehfest
// ) {
//
//   time = std::log(2.0) / time.array();
//
//   Eigen::VectorXd v = stehfest_v(n_stehfest);
//   Eigen::MatrixXd p = stehfest_p(time, n_stehfest);
//
//   std::vector<double> p_vec(p.data(), p.data() + p.size());
//
//   for (auto &out : p_vec)
//     out = theis_large_diameter_laplace(out, radius, radius_well, t_1, s_1);
//
//   p = Eigen::Map<Eigen::MatrixXd>(p_vec.data(), p.rows(), p.cols());
//
//   p.array() *= flow_rate;
//
//   Eigen::VectorXd ret = (p * v).array() * time.array();
//
//   return Rcpp::List::create(Rcpp::Named("theis_large_diameter") = ret);
//
// }
//
//
//
//
// /*** R
//
// x <- sort(abs(rnorm(1000)))
// y <- rnorm(10)
// x2 <- log(2.0) / x
//
// hydrorecipes:::stehfest_v(12L)
// hydrorecipes:::stehfest_p(as.numeric(1:n), 12L)
//
// barker_herbert(as.numeric(1:n), 100.0, 200.0, 1.0, 1e-3, 1e-3, 1e-5, 1e-5)
//
// (hydrorecipes:::stehfest_barker_herbert(as.numeric(1:n),1.0, 100.0, 200.0, 1.0, 1e-3, 1e-3, 1e-5, 1e-5, 12L))
//
//
// n <- 1000000
// bench::mark(
// (hydrorecipes:::stehfest_barker_herbert(as.numeric(1:n), 1.0,
//                                     100.0, 200.0,  1e-3, 1e-3, 1e-5, 1e-5, 12L)[[1]])
// )
//
//
// n <- 100000
// bench::mark(
//   x <- (hydrorecipes:::stehfest_theis_large_diameter(
//     as.numeric(1:n)/1,
//     0.1,
//     100.0,
//     0.2,
//     1e-4,
//     1e-5,
//     12L)[[1]])
// )
// plot(x, type = 'l', log = 'xy')
// points(hydrorecipes:::grf_time(radius = 100.0,
//                            1e-5,
//                            1e-4,
//                            1.0,
//                            as.numeric(1:n)/1,
//                            flow_rate = rep(.1, n) ,
//                            2.0)[[1]], type = 'l',
//        col = 'red')
// abline(a = 0, b = 1)
// hydrorecipes:::grf_time()
// n <- 100
// bench::mark(
// barker_herbert(as.numeric(1:n), 100.0, 200.0, 1.0, 1e-3, 1e-3, 1e-5, 1e-5),
// (hydrorecipes:::stehfest_barker_herbert(as.numeric(1:n),
//                                     1.0,
//                                     100.0, 200.0,
//                                     1e-3, 1e-3, 1e-5, 1e-5, 12L)[[1]])
// )
// */
