#include "frecipes.h"

using namespace boost::math;

struct ParallelFracturesSolute {
  Eigen::VectorXd tau;
  double c_0;
  double lambda;
  double nu;
  double z;
  double x;
  double a;
  double g;
  double b;
  double B;
  double sigma;
  double k2;

  ParallelFracturesSolute(Eigen::VectorXd time,
                          double z,
                          double x,
                          double c_0,
                          double b,
                          double B,
                          double v,
                          double t_half,
                          double D_star,
                          double alpha_l,
                          double tortuosity, // τ
                          double k_m,
                          double k_f,
                          double rho_b,
                          double theta) :
    z(z), x(x), c_0(c_0), B(B) {

    // laplace variable
    tau = std::log(2.0) / time.array();

    // D'
    double D_prime = tortuosity * D_star;

    // D
    double D = alpha_l * v + D_star;

    // λ
    // equation 3
    lambda = std::log(2.0) / t_half;

    // R = 1 + K_f / b
    // equation 4
    double R = 1.0 + k_f / b;

    // R' = 1 + ρ_b / θ K_m
    double R_prime = 1.0;
    if (theta > 0.0) {
      R_prime = 1.0 + (rho_b / theta) * k_m;
    }

    // equation 13
    g = sqrt(R_prime / D_prime);

    // equation 14
    sigma = g * (0.5 * (B - b));

    // equation 16
    a = ((0.5 * b) * R) / (theta * sqrt(R_prime * D_prime));

    // equation 17
    nu = v / (2.0 * D);
    k2 = 4.0 * R * D / (v * v);

  };

  double lp(double p) {
    // Rcout << "The value is p " << p << std::endl;
    p = p + lambda;

    double sqrt_p = sqrt(p);

    // equation 19
    double c_fracture =
      c_0 / (p - lambda) *
      exp(nu * z) *
      (exp(-(nu * z) * sqrt(1.0 + k2 * (p + sqrt_p / a * std::tanh(sigma * sqrt_p)))));

    if(x <= (0.5 * b)) {
      return(c_fracture);
    }

    return(
      c_fracture *
        (std::cosh(g * sqrt_p * ((0.5 * B) - x)) /
         std::cosh(sigma * sqrt_p)));
  };
  Eigen::VectorXd lp_2(double p, Eigen::VectorXd zs, Eigen::VectorXd xs) {

    p = p + lambda;

    double sqrt_p = sqrt(p);
    double p_term = sqrt(1.0 + k2 * (p + sqrt_p / a * std::tanh(sigma * sqrt_p)))
    double z_term, c_fracture;

    int n_z zs.size();
    int n_x xs.size();
    Eigen::VectorXd c_fracture(n_z);
    Eigen::RowVectorXd x_term(n_x);

    // Each z location
    for(int i = 0; i < n_z; ++i) {
      z_term = exp(nu * zs[i]);
      c_fracture[i] = c_0 / (p - lambda) * z_term * (1.0 / z_term * term));
    }

    // Each x location
    for(int i = 0; i < n_x; ++i) {
      x_term[i] = (std::cosh(g * sqrt_p * ((0.5 * B) - x)) /
                   std::cosh(sigma * sqrt_p)));
    }

    Eigen::MatrixXd ret_mat = c_fracture * x_term;
    VectorXd ret(Map<VectorXd>(ret_mat.data(), ret_mat.cols()*A.rows()));

    return(ret);

  };


};

// for each p matrix



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

// [[Rcpp::export]]
Rcpp::List parallel_fractures_chem(
    Eigen::VectorXd time,
    double z,
    double x,
    double c_0,
    double b,
    double B,
    double v,
    double t_half,
    double D_star,
    double alpha_l,
    double tortuosity, // τ
    double k_m,
    double k_f,
    double rho_b,
    double theta,
    int n_stehfest) {

  ParallelFracturesSolute fracture(time,
                                   z,
                                   x,
                                   c_0,
                                   b,
                                   B,
                                   v,
                                   t_half,
                                   D_star,
                                   alpha_l,
                                   tortuosity, // τ
                                   k_m,
                                   k_f,
                                   rho_b,
                                   theta);


    return(stehfest(fracture, n_stehfest));

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

# bench::mark(GCD(as.integer(times)))
flow_rate <- rep(Q, n)

bench::mark(
  a <- frecipes:::hantush_jacob_laplace(times, cc, rw, Tr, S, Q, prec, n_terms),
  b <- frecipes:::hantush_jacob(times, flow_rate, rw,S,Tr, lab, prec),
  # c <- frecipes:::hantush_jacob_quad(times, lab, rw, Tr, S, Q, 1e-16),
  # c <- frecipes:::barker_herbert(times, c, rw, Tr, S, Q, prec, 12L),
  check = FALSE
)



# ParallelFracturesSolute fracture(time,
#                                  z,
#                                  x,
#                                  c_0,
#                                  b,
#                                  B,
#                                  v,
#                                  t_half,
#                                  D_star,
#                                  alpha_l,
#                                  tortuosity, // τ
#                                  k_m,
#                                  k_f,
#                                  rho_b,
#                                  theta);

lap <- frecipes:::parallel_fractures_chem(
  seq(1, 20000, length.out = 1000), # time
  5.86/4, # z
  0.0, # x
  1.0, # c_0
  100 * 1e-6, # aperture
  0.5, # spacing
  0.0075, # v
  365 * 12.35, # t_half
  1.6e-9*86400, # D_star
  0.1, # alpha_l
  0.1, # tortuosity
  0.0, # k_m
  0.0, # k_f
  0.0, # rho_b
  0.01, # theta porosity
  18L)[[1]]


plot(lap, type = "l", ylim = c(0,1))


lap <- frecipes:::parallel_fractures_chem(
  seq(1, 10000, length.out = 10000), # time
  0.5, # z
  0.1, # x
  1.0, # c_0
  100 * 1e-6, # aperture
  0.5, # spacing
  0.1, # v
  365 * 12.35, # t_half
  1.6e-9*86400, # D_star
  0.1, # alpha_l
  0.1, # tortuosity
  0.0, # k_m
  0.0, # k_f
  0.0, # rho_b
  0.1, # theta porosity
  18L)[[1]]

points(lap, type = "l", col = 'red')

lap

rangelaprange(abs(a[[1]]-c))
range(abs(-b[[1]]-c))


plot(c, type = 'l')
points(-b[[1]], type = 'l', col = 'blue')
points(a[[1]], type = 'l', col = 'red')


e1 <- exp(2.3)
bench::mark(
  exp(-2.3),
  1/exp(2.3),
  1/e1
)
*/
