#include "hydrorecipes.h"


// [[Rcpp::export]]
Eigen::MatrixXd gwr_p(
    Eigen::VectorXd time,
    unsigned int n_gwr
) {

  Eigen::RowVectorXd v = Eigen::RowVectorXd::LinSpaced(n_gwr * 2, 1.0, (double)n_gwr * 2.0);

  return(time * v);

}




// [[Rcpp::export]]
double barker_herbert_impulse(double p,
                              double radius,
                              double radius_patch,
                              double t_1,
                              double t_2,
                              double s_1,
                              double s_2) {

  double N = sqrt(s_1 * p / t_1);
  double A = sqrt(s_2 * p / t_2);

  double ct = (t_2 / t_1) * (A / N);

  double bi_n0 = std::cyl_bessel_i(0.0, N * radius_patch);
  double bk_a0 = std::cyl_bessel_k(0.0, A * radius_patch);
  double bk_n0 = std::cyl_bessel_k(0.0, N * radius_patch);
  double bi_n1 = std::cyl_bessel_i(1.0, N * radius_patch);
  double bk_a1 = std::cyl_bessel_k(1.0, A * radius_patch);
  double bk_n1 = std::cyl_bessel_k(1.0, N * radius_patch);

  double denom = (ct * bi_n0 * bk_a1 + bi_n1 * bk_a0) * p;

  double term_1 = (bk_n1 * bk_a0  - bk_a1 * bk_n0 * ct);
  double term_2 = (bk_n0 * bi_n1 + bk_n1 * bi_n0);


  double drawdown;

  // if radius is inside the patch
  if (radius <= radius_patch) {

    drawdown = std::cyl_bessel_k(0.0, N * radius) / p  +
     (term_1 * std::cyl_bessel_i(0.0, N * radius)) / denom;

    return(drawdown);

  }

  drawdown = term_2 * std::cyl_bessel_k(0.0, A * radius) / denom;

  return(drawdown);

}



// [[Rcpp::export]]
Rcpp::List gwr_barker_herbert(
     Eigen::VectorXd time,
     double flow_rate,
     double radius,
     double radius_patch,
     double t_1,
     double t_2,
     double s_1,
     double s_2,
     unsigned int n_gwr
) {

  Eigen::VectorXd ret(time.size());
  bool broken = false;
  // time to tau
  time = std::log(2.0) / time.array();
  Eigen::MatrixXd p = gwr_p(time, n_gwr);
  std::vector<double> p_vec(p.data(), p.data() + p.size());
  double sm;


  for (auto &out : p_vec)
    out = barker_herbert_impulse(out, radius, radius_patch, t_1, t_2, s_1, s_2);

  p = Eigen::Map<Eigen::MatrixXd>(p_vec.data(), p.rows(), p.cols());
  p.array() *= flow_rate / (2.0 * M_PI * t_1);

  unsigned int m1 = n_gwr;


  Eigen::VectorXd g0 = Eigen::VectorXd::Zero(m1 + 1);
  for (unsigned int kk = 0; kk < time.size(); ++kk) {
    double tau = time[kk];
    for (unsigned int n = 1; n <= n_gwr; ++n) {
      sm = 0.0;
      for (unsigned int i = 0; i <= n; ++i) {
        sm = R::choose(n, i) * pow(-1.0, i) * p(kk, n+i);  // can we do this convolution faster with fft (overlap add/save)? or matrix multiply?
      }
      g0[n] = tau * std::tgamma(2 * n + 1) /
        (std::tgamma(n + 1) * std::tgamma(n)) * sm;
    }
  Eigen::VectorXd gm = Eigen::VectorXd::Zero(m1 + 1);

  Eigen::VectorXd gp = Eigen::VectorXd::Zero(m1 + 1);
  double best = g0[m1];
  double expr;
  for (unsigned int k = 0; k < m1 - 2; ++k) {
    for (unsigned int n = m1 -2 -k; n > 0; --n) {
      expr = g0[n + 2] - g0[n + 1];
      if (expr == 0){
        broken = true;
        break;
      }
      expr = gm[n + 2] + (k + 1) / expr;
      gp[n + 1] = expr;
      if ((k / 2 == 1) && (n == (m1 - 2 - k))){
        best = expr;
      }
    }
    if (broken) {
      break;
    }
    for (unsigned int n = 0; n > m1-k; ++n) {
      gm[n + 1] = g0[n + 1];
      g0[n + 1] = gp[n + 1];
    }
  }
  ret[kk] = best;
  }

  // Eigen::VectorXd ret = (p * v).array() * time.array();

  return Rcpp::List::create(Rcpp::Named("barker_herbert") = ret);

}

/*** R
n <- 100
bench::mark(

(hydrorecipes:::gwr_barker_herbert(as.numeric(1:n),
                                    1.0,
                                    100.0, 200.0,
                                    1e-3, 1e-3, 1e-5, 1e-5, 12L)[[1]])
)

# (hydrorecipes:::stehfest_barker_herbert(as.numeric(1:n),
#                                     1.0,
#                                     100.0, 200.0,
#                                     1e-3, 1e-3, 1e-5, 1e-5, 12L)[[1]])

*/
