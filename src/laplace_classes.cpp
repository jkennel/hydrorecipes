#include "frecipes.h"

using namespace boost::math;

struct HantushJacob {
  Eigen::VectorXd tau;
  double c;
  double r;
  double Tr;
  double S;
  double Q;
  double prec;
  HantushJacob(Eigen::VectorXd time,
               double c,
               double r,
               double Tr,
               double S,
               double Q,
               double prec) : c(c), r(r), Tr(Tr), S(S), Q(Q), prec(prec) {
    tau = std::log(2.0) / time.array();
  };
  double lp(double p) {
    double w = (S * c * p + 1.0) / (c * Tr);
    return (-Q / (2.0 * M_PI * Tr * p) * cyl_bessel_k(0.0, r * sqrt(w)));
  };
};

struct BarkerHerbert {
  Eigen::VectorXd tau;
  double radius;
  double radius_patch;
  double t_1;
  double t_2;
  double s_1;
  double s_2;
  double Q;
  double prec;
  BarkerHerbert(Eigen::VectorXd time,
                double radius,
                double radius_patch,
                double t_1,
                double t_2,
                double s_1,
                double s_2,
                double Q,
                double prec) : radius(radius), radius_patch(radius_patch),
                t_1(t_1), t_2(t_2), s_1(s_1), s_2(s_2), Q(Q), prec(prec) {
    tau = std::log(2.0) / time.array();
  };
  double lp(double p) {
    double N = sqrt(s_1 * p / t_1);
    double A = sqrt(s_2 * p / t_2);

    double ct = (t_2 / t_1) * (A / N);

    double bi_n0 = cyl_bessel_i(0.0, N * radius_patch);
    double bk_a0 = cyl_bessel_k(0.0, A * radius_patch);
    double bk_n0 = cyl_bessel_k(0.0, N * radius_patch);
    double bi_n1 = cyl_bessel_i(1.0, N * radius_patch);
    double bk_a1 = cyl_bessel_k(1.0, A * radius_patch);
    double bk_n1 = cyl_bessel_k(1.0, N * radius_patch);

    double denom = (ct * bi_n0 * bk_a1 + bi_n1 * bk_a0) * p;

    double term_1 = (bk_n1 * bk_a0  - bk_a1 * bk_n0 * ct);
    double term_2 = (bk_n0 * bi_n1 + bk_n1 * bi_n0);

    double drawdown;

    // if radius is inside the patch
    if (radius <= radius_patch) {

      drawdown = cyl_bessel_k(0.0, N * radius) / p  +
        (term_1 * cyl_bessel_i(0.0, N * radius)) / denom;

      return(drawdown);

    }

    // if radius is outside the patch
    drawdown = term_2 * cyl_bessel_k(0.0, A * radius) / denom;

    return(drawdown * Q / (2.0 * M_PI * t_1));

  }
};

// [[Rcpp::export]]
Eigen::RowVectorXd stehfest_v(
    int n
) {

  assert(n < 20);
  assert(n % 2 == 0);

  int n_div_2 = n / 2;
  int s, e;
  double z;

  Eigen::RowVectorXd fact(n+1);
  Eigen::RowVectorXd v(n);

  // calculate factorials
  for (unsigned int i = 0; i <= n; ++i) {
    fact[i] = factorial<double>(i);
  }

  for (int i = 0; i < n; ++i) {
    if (i > n_div_2 - 1) {
      s = n_div_2 - 1;
    } else {
      s = i;
    }
    e = (int)(i / 2);

    z = 0.0;
    for (int k = s; k >= e; --k) {
      z += (pow((double)k + 1.0, n_div_2) * fact[2 * k + 2]) /
        (fact[n_div_2 - k - 1] * fact[k + 1] * fact[k] *
          fact[i - k] * fact[2 * k - i + 1]);

    }

    v[i] = pow(-1.0, n_div_2 + i + 1.0) * z;

  }

  return(v);

}


// [[Rcpp::export]]
Eigen::MatrixXd stehfest_p(
    Eigen::VectorXd time,
    int n_stehfest
) {

  Eigen::RowVectorXd v = Eigen::RowVectorXd::LinSpaced(n_stehfest, 1.0, (double)n_stehfest);

  return(time * v);

}


template <typename T>
Rcpp::List stehfest(T& well, int n_stehfest) {

  Eigen::VectorXd v = stehfest_v(n_stehfest);
  Eigen::MatrixXd p = stehfest_p(well.tau, n_stehfest);

  std::vector<double> p_vec(p.data(), p.data() + p.size());

  for (auto &out : p_vec)
    out = well.lp(out);

  p = Eigen::Map<Eigen::MatrixXd>(p_vec.data(), p.rows(), p.cols());
  Eigen::VectorXd ret = (p * v).array() * (well.tau).array();

  return Rcpp::List::create(Rcpp::Named("lap") = ret);
}

// [[Rcpp::export]]
Rcpp::List hantush_jacob_laplace(
    Eigen::VectorXd time,
    double c,
    double r,
    double Tr,
    double S,
    double Q,
    double prec,
    int n_stehfest) {

  HantushJacob well(time, c, r, Tr, S, Q, prec);

  return(stehfest(well, n_stehfest));

}

// [[Rcpp::export]]
Rcpp::List barker_herbert(
    Eigen::VectorXd time,
    double radius,
    double radius_patch,
    double t_1,
    double t_2,
    double s_1,
    double s_2,
    double Q,
    double prec,
    int n_stehfest) {

  BarkerHerbert well(time, radius, radius_patch,
                     t_1, t_2, s_1, s_2, Q, prec);

  return(stehfest(well, n_stehfest));

}



/*** R
n <- 10000

Tr = 200 # transmissivity of aquifer, m^2/d
S = 0.0005 # storage coefficient of aquifer, -
cc = 1000 # resistance of leaky layer, d
Q = 800 # discharge of well, m^3/d
rw = 2500 # radius of well, m
lab = sqrt(cc * Tr)
n_terms <- 8L
prec = 1e-5
times <- round(2*c(seq(0, 3000000, by = 60), 60), 0)

bench::mark(GCD(as.integer(times)))
flow_rate <- rep(Q, n)

bench::mark(
  a <- frecipes:::hantush_jacob_laplace(times, cc, rw, Tr, S, Q, prec, n_terms),
  b <- frecipes:::hantush_jacob(times, flow_rate, rw,S,Tr, lab, prec),
  c <- frecipes:::hantush_jacob_quad(times, lab, rw, Tr, S, Q, 1e-16),
  # c <- frecipes:::barker_herbert(times, c, rw, Tr, S, Q, prec, 12L),
  check = FALSE
)


range(abs(a[[1]]-c))
range(abs(-b[[1]]-c))


plot(c, type = 'l')
points(-b[[1]], type = 'l', col = 'blue')
points(a[[1]], type = 'l', col = 'red')

*/
