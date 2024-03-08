//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
//
// Laplace Domain methods  -----------------------------------------------------
//
//
// These codes have not been optimize for large datasets and I'm not sure of the
// current value.
//
//  - Some results may oscillate
//  - Should arbitrary precision be used?
//  - Calculate with log time spacing followed by interpolation to speed up?
//
//
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#include "frecipes.h"


//' @title
//' bessel_k_cplx
//'
//' @description
//' Modified Bessel function of first kind order 1
//'
//' @param x \code{numeric} value to evaluate
//' @param nu \code{numeric} value to evaluate
//' @param expon_scaled \code{boolean} value to evaluate
//' @param n_seq \code{nseq} value to evaluate
//'
//' @return bessel function result
//'
//'
//' @export
// [[Rcpp::export]]
Rcpp::ComplexMatrix bessel_k_cplx(const Rcpp::ComplexMatrix &x,
                                  double nu,
                                  bool expon_scaled,
                                  size_t n_seq)
{

  unsigned int n = x.size();
  // Obtain environment containing function
  Rcpp::Environment Bessel("package:Bessel");

  // Make function callable from C++
  Rcpp::Function bessel_k = Bessel["BesselK"];

  if(n_seq == 1) {
    Rcpp::ComplexVector v = bessel_k(x, nu, expon_scaled, n_seq);
    Rcpp::ComplexMatrix m(n, 1, v.begin());
    return(m);
  }

  return(bessel_k(x, nu, expon_scaled, n_seq));

}


struct CooperBredehoeftPapadopulos
{
  Eigen::VectorXd tau;
  double r;
  double r_c;
  double r_w;
  double Tr;
  double S;
  double h_0;
  double alpha;
  CooperBredehoeftPapadopulos(
    Eigen::VectorXd time,
    double r,
    double r_c,
    double r_w,
    double Tr,
    double S,
    double h_0) : r(r), r_w(r_w), r_c(r_c), Tr(Tr), S(S), h_0(h_0)
  {
    alpha = (r_w * r_w * S) / (r_c * r_c);
    tau = std::log(2.0) / time.array();
  };
  double lp(double p)
  {
    // Rcpp::Rcout << "The value is p " << p << std::endl;

    if (std::isinf(p)) {
      return(p);
    }

    double q = sqrt((p * S) / Tr);

    double bk_rw = std::cyl_bessel_k(0.0, r_w * q);
    double bk_r = bk_rw;
    double bk_rw_1 = std::cyl_bessel_k(1.0, r_w * q);

    if (r > r_w) {
      bk_r = std::cyl_bessel_k(0.0, r * q);
    }

    double f_p = (r_c * S * h_0 * bk_r) /
      ((Tr * q) * ((r_w * q * bk_rw) + (2.0 * alpha * bk_rw_1)));

    return(f_p);

  };
};



struct ParallelFracturesHeat
{
  Eigen::VectorXd time;
  Eigen::VectorXd gamma; // gamma
  Eigen::VectorXd tau;
  Eigen::VectorXd z;
  Eigen::VectorXd x;
  Eigen::VectorXd del_c;
  Eigen::VectorXd mean_time;
  double t_0;
  double lambda;
  double nu;
  double a;
  double g;
  double b;
  double B;
  double sigma;
  double k2;

  ParallelFracturesHeat(Eigen::VectorXd time,
                        Eigen::VectorXd z,
                        Eigen::VectorXd x,
                        Eigen::VectorXd temperature_influent,
                        Eigen::VectorXd time_influent,
                        double t_0,
                        double b,
                        double B,
                        double v,
                        double lambda_fracture, // thermal conductivity of water
                        double lambda_matrix,   // thermal conductivity of solids
                        double spec_heat_w,     // specific heat of water
                        double spec_heat_s,     // specific heat of solids
                        double rho_w,           // density of water
                        double rho_s,           // density of solids
                        double theta,           // porosity
                        unsigned int n_terms) : time(time),
                        z(z), x(x),
                        t_0(t_0),
                        B(B), b(b)

  {

    lambda = 0.0;

    // laplace variable Stehfest
    tau = std::log(2.0) / time.array();

    // laplace variable Cohen
    gamma = (2.0 / 3.0) *
      (std::floor(n_terms / 1.31) + log(10.0) + (2.0 * time.array()).log());

    unsigned int n_c = temperature_influent.size();
    unsigned int n_t = time_influent.size();

    if (n_c != n_t)
    {
      Rcpp::stop("number of influent times and influent concentrations must be equal");
    }

    del_c = temperature_influent.array();

    mean_time = time_influent;

    if (n_c > 1)
    {
      for (unsigned int i = 1; i < n_c; ++i)
      {
        mean_time[i] = (time_influent[i] + time_influent(i - 1)) / 2.0;
        del_c[i] = (temperature_influent[i]) - (temperature_influent(i - 1));
      }
    }

    // D'
    double D_prime = lambda_matrix /
      ((theta * spec_heat_w * rho_w) + ((1.0 - theta) * spec_heat_s * rho_s));

    // D
    double D = lambda_fracture / (spec_heat_w  * rho_w);

    double R = 1.0;
    double R_prime = 1.0;

    g = sqrt(R_prime / D_prime);
    sigma = g * (0.5 * (B - b));

    // equation 16
    a = (lambda_matrix * g) / (spec_heat_w * rho_w * b * 0.5);

    // equation 17
    nu = v / (2.0 * D);
    k2 = 4.0 * R * D / (v * v);
  };

  double cohen_xy(std::complex<double> p, double t_p, double z, double x)
  {
    // Inflow
    std::complex<double> sum(0.0, 0.0);
    std::complex<double> term(0.0, 0.0);
    // Rcout << "The value is p " << p << std::endl;

    for (unsigned int i = 0; i < del_c.size(); ++i)
    {
      if (mean_time[i] < t_p)
      {
        if (mean_time[i] <= 0.0)
        {
          term = del_c[i] / (p);
        }
        else
        {
          term = (del_c[i] / p) * exp(-p * mean_time[i]);
        }
      }
      else
      {
        term = std::complex<double>(0.0, 0.0);
      }
      sum += term;
    }

    double n_cut = 710.0;

    std::complex<double> sqrt_p = std::sqrt(p);
    std::complex<double> t_0_p = t_0 / p;
    std::complex<double> sigma_sqrt_p = sigma * sqrt_p;
    std::complex<double> tan_h;
    std::complex<double> xx = (g * sqrt_p * ((0.5 * B) - x));
    std::complex<double> xterm(1.0, 0.0);
    std::complex<double> diff = xx - sigma_sqrt_p;

    if(sigma_sqrt_p.real() > n_cut) {
      tan_h = std::complex<double>(1.0, 0.0);
      if(diff.real() > -n_cut & diff.real() < n_cut & x > 0.5 * b) {
        xterm = std::exp(diff);
      } else if (x > 0.5 * b) {
        xterm = std::complex<double>(0.0, 0.0);
      }
    } else if (sigma_sqrt_p.real() < -n_cut) {
      tan_h = std::complex<double>(-1.0, 0.0);
      if(diff.real() > -n_cut & diff.real() < n_cut & x > 0.5 * b) {
        xterm = std::exp(diff);
      } else if (x > 0.5 * b) {
        xterm = std::complex<double>(0.0, 0.0);
      }
    } else {
      tan_h = std::tanh(sigma_sqrt_p);
      // Rcout << "The value is tanh " << tan_h << std::endl;

      if (x > 0.5 * b) {
        xterm = std::cosh(xx) / std::cosh(sigma_sqrt_p);
      }
    }

    std::complex<double> root = nu *
      (1.0 - std::sqrt(1.0 + k2 * (p + sqrt_p * a * tan_h)));

    std::complex<double> c_fracture = (sum - t_0_p) * std::exp(root * z) + t_0_p;

    return ((c_fracture - t_0_p) * xterm + t_0_p).real();
  };
};

struct ParallelFracturesSolute
{
  Eigen::VectorXd time;
  Eigen::VectorXd gamma; // gamma
  Eigen::VectorXd tau;
  Eigen::VectorXd z;
  Eigen::VectorXd x;
  Eigen::VectorXd del_c;
  Eigen::VectorXd mean_time;
  double c_0;
  double lambda;
  double nu;
  double a;
  double g;
  double b;
  double B;
  double sigma;
  double k2;

  ParallelFracturesSolute(Eigen::VectorXd time,
                          Eigen::VectorXd z,
                          Eigen::VectorXd x,
                          Eigen::VectorXd concentration_influent,
                          Eigen::VectorXd time_influent,
                          double c_0,
                          double b,
                          double B,
                          double v,
                          double alpha_l,
                          double D_star,
                          double k_f,
                          double k_m,
                          double t_half,
                          double rho_b,
                          double theta,
                          double tortuosity, // τ
                          unsigned int n_terms) : time(time),
                          z(z), x(x),
                          c_0(c_0),
                          B(B), b(b)

  {

    // laplace variable Stehfest
    tau = std::log(2.0) / time.array();

    // laplace variable Cohen
    gamma = (2.0 / 3.0) *
      (std::floor(n_terms / 1.31) + log(10.0) + (2.0 * time.array()).log());

    unsigned int n_c = concentration_influent.size();
    unsigned int n_t = time_influent.size();

    if (n_c != n_t)
    {
      Rcpp::stop("number of influent times and influent concentrations must be equal");
    }

    del_c = concentration_influent;
    mean_time = time_influent;

    if (n_c > 1)
    {
      for (unsigned int i = 1; i < n_c; ++i)
      {
        mean_time[i] = (time_influent[i] + time_influent[i - 1]) / 2.0;
        del_c[i] = concentration_influent[i] - concentration_influent[i - 1];
      }
    }

    // D'
    double D_prime = tortuosity * D_star;

    // D
    double D = alpha_l * v + D_star;

    // λ
    // equation 3
    lambda = std::log(2.0) / t_half;

    // R = 1 + K_f / b
    // equation 4
    double R = 1.0 + k_f / (b / 2);

    // R' = 1 + ρ_b / θ K_m
    double R_prime = 1.0;
    if (theta > 0.0)
    {
      R_prime = 1.0 + (rho_b / theta) * k_m;
    }

    // equation 13
    g = sqrt(R_prime / D_prime);

    // equation 14
    sigma = g * (0.5 * (B - b));

    // equation 16
    a = (theta * sqrt(R_prime * D_prime)) / ((0.5 * b) * R);

    // equation 17
    nu = v / (2.0 * D);
    k2 = 4.0 * R * D / (v * v);
  };

  double lp_xy(double p, double t_p, double z, double x)
  {
    // Inflow
    double sum = 0;
    double term = 0;

    for (unsigned int i = 0; i < del_c.size(); ++i)
    {

      if (mean_time[i] < t_p)

      {
        if (mean_time[i] <= 0.0)
        {
          term = del_c[i] / (p);
        }
        else
        {
          term = (del_c[i] / p) * exp(-p * mean_time[i]);
        }
      }
      else
      {
        term = 0.0;
      }
      sum += term;
    }

    // Rcout << "The value is p " << p << std::endl;
    p = p + lambda;

    if (x > 0.5 * B)
    {
      Rcpp::stop("x should be less than 0.5 B");
    }
    double n_cut = 710.0;
    double sqrt_p = sqrt(p);
    double sigma_sqrt_p = sigma * sqrt_p;
    double c_0_p = c_0 / p;

    double tan_h;
    double xx = (g * sqrt_p * ((0.5 * B) - x));
    double xterm = 1.0;
    double diff = xx - sigma_sqrt_p;

    // check bounds for tanh and cosh
    if(sigma_sqrt_p > n_cut) {
      tan_h = 1.0;
      if(diff > -n_cut & diff < n_cut) {
        xterm = std::exp(diff);
      } else {
        xterm = 0.0;
      }
    } else if (sigma_sqrt_p < -n_cut) {
      tan_h = -1.0;
      if(diff > -n_cut & diff < n_cut) {
        xterm = std::exp(diff);
      } else {
        xterm = 0.0;
      }
    } else {
      tan_h = std::tanh(sigma_sqrt_p);
      xterm = std::cosh(xx) / std::cosh(sigma_sqrt_p);
    }

    double root = nu * (1.0 - sqrt(1.0 + k2 * (p + sqrt_p * a * tan_h)));
    double c_fracture = (sum - c_0_p) * std::exp(root * z) + c_0_p;

    return ((c_fracture - c_0_p) * xterm + c_0_p);
  };
  double cohen_xy(std::complex<double> p, double t_p, double z, double x)
  {
    // Inflow
    std::complex<double> sum(0.0, 0.0);
    std::complex<double> term(0.0, 0.0);

    // Rcpp::Rcout << "==================" << t_p << std::endl;
    for (unsigned int i = 0; i < del_c.size(); ++i)
    {
      // Rcpp::Rcout << "The value is mean_time" << mean_time[i] << std::endl;

      if (mean_time[i] < t_p)

      {
        if (mean_time[i] <= 0.0)
        {
          term = del_c[i] / (p);
        }
        else
        {
          term = (del_c[i] / p) * exp(-p * mean_time[i]);
        }
      }
      else
      {
        term = std::complex<double>(0.0, 0.0);
      }
      sum += term;
    }

    // Rcpp::Rcout << "The value is p " << p << std::endl;
    // Rcpp::Rcout << "The value is sum " << sum << std::endl;
    // Rcpp::Rcout << "The value is term " << term << std::endl;

    // Rcout << "The value is p " << p << std::endl;
    p = p + lambda;
    if (x > 0.5 * B)
    {
      Rcpp::stop("x should be less than 0.5 B");
    }
    double n_cut = 710.0;

    std::complex<double> sqrt_p = std::sqrt(p);
    std::complex<double> c_0_p = c_0 / p;
    std::complex<double> sigma_sqrt_p = sigma * sqrt_p;
    std::complex<double> tan_h;
    std::complex<double> xx = (g * sqrt_p * ((0.5 * B) - x));
    std::complex<double> xterm(1.0, 0.0);
    std::complex<double> diff = xx - sigma_sqrt_p;

    if(sigma_sqrt_p.real() > n_cut) {
      tan_h = std::complex<double>(1.0, 0.0);
      if(diff.real() > -n_cut & diff.real() < n_cut & x > 0.5 * b) {
        xterm = std::exp(diff);
      } else if (x > 0.5 * b) {
        xterm = std::complex<double>(0.0, 0.0);
      }
    } else if (sigma_sqrt_p.real() < -n_cut) {
      tan_h = std::complex<double>(-1.0, 0.0);
      if(diff.real() > -n_cut & diff.real() < n_cut & x > 0.5 * b) {
        xterm = std::exp(diff);
      } else if (x > 0.5 * b) {
        xterm = std::complex<double>(0.0, 0.0);
      }
    } else {
      tan_h = std::tanh(sigma_sqrt_p);
      if (x > 0.5 * b) {
        xterm = std::cosh(xx) / std::cosh(sigma_sqrt_p);
      }
    }


    std::complex<double> root = nu *
      (1.0 - std::sqrt(1.0 + k2 * (p + sqrt_p * a * tan_h)));
    std::complex<double> c_fracture = (sum - c_0_p) * std::exp(root * z) + c_0_p;

    // Rcpp::Rcout << "-------------" << nu << std::endl;
    // Rcpp::Rcout << "The value is sigma_sqrt_p " << sigma_sqrt_p << std::endl;
    // Rcpp::Rcout << "The value is nu " << nu << std::endl;
    // Rcpp::Rcout << "The value is a " << a << std::endl;
    // Rcpp::Rcout << "The value is sigma " << sigma << std::endl;
    // Rcpp::Rcout << "The value is k2 " << k2 << std::endl;
    // Rcpp::Rcout << "The value is sqrt_p " << sqrt_p << std::endl;
    // Rcpp::Rcout << "The value is c_0_p " << c_0_p << std::endl;
    // Rcpp::Rcout << "The value is root " << root << std::endl;
    // Rcpp::Rcout << "The value is sqrt_p " << sqrt_p << std::endl;
    // Rcpp::Rcout << "The value is xterm " << sqrt_p << std::endl;
    // Rcpp::Rcout << "The value is c_fracture " << sqrt_p << std::endl;



    // if (x > 0.5 * b) // in the matrix
    // {
    //   if((g * sqrt_p * ((0.5 * B) - x)).real() < n_cut & sigma_sqrt_p.real() < n_cut) {
    //   } else {
    //     xterm = std::complex<double>(0.0, 0.0);
    //   }
    // }

    return ((c_fracture - c_0_p) * xterm + c_0_p).real();
  };
};



struct JacobLohman
{
  Eigen::VectorXd tau;
  double s;
  double r;
  double Tr;
  double S;
  double prec;
  JacobLohman(Eigen::VectorXd time,
               double s,
               double r,
               double Tr,
               double S,
               double prec) : r(r), Tr(Tr), s(s), S(S), prec(prec)
  {
    tau = std::log(2.0) / time.array();

  };
  double lp(double p)
  {
    double w = sqrt(p * S / Tr);
    double dbar = s * 1.0 / p;
    return (dbar * w * std::cyl_bessel_k(1.0, r * w) / std::cyl_bessel_k(0.0, r * w));
  };
  // double cohen(std::complex<double> p)
  // {
  //   std::complex<double> w = sqrt(p * S / Tr);
  //   std::complex<double> dbar = s * 1.0 / p;
  //   return (dbar * w * std::cyl_bessel_k(1.0, r * w) / std::cyl_bessel_k(0.0, r * w)).real();
  // };
};

struct PapadopulosCooper
{
  Eigen::VectorXd tau;
  double Q;
  double r;
  double r_c;
  double r_w;
  double Tr;
  double S;
  double fact;
  double prec;
  PapadopulosCooper(
    Eigen::VectorXd time,
    double Q,
    double r,
    double r_c,
    double r_w,
    double Tr,
    double S,
    double prec) : Q(Q), r(r), r_w(r_w), r_c(r_c), Tr(Tr), S(S), prec(prec)
  {
    fact = 2.0 * M_PI * Tr;
    tau = std::log(2.0) / time.array();
  };
  double lp(double p)
  {
    double w = sqrt(p * S / Tr);
    double dbar = Q * 1.0 / p;
    double a1 = r_w * w;
    double a2 = r * w;

    double f_p = Q * std::cyl_bessel_k(0, a2) /
      (M_PI * p * ((r_c * r_c * p * std::cyl_bessel_k(0.0, a1)) +
      (2 * r_w * Tr * w * std::cyl_bessel_k(1.0, a1))));

    return(f_p);

    // double term_1 = (r_c * r_c * p / (2.0 * Tr)) *
    //   std::cyl_bessel_k(0,a1) / std::cyl_bessel_k(0, a2);
    // double term_2 = r_w * w * std::cyl_bessel_k(1.0, a1)/std::cyl_bessel_k(0, a2);
    // return (dbar / (term_1 + term_2));

  };
};



struct HantushJacob
{
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
               double prec) : c(c), r(r), Tr(Tr), S(S), Q(Q), prec(prec)
  {

    tau = std::log(2.0) / time.array();
  };
  double lp(double p)
  {
    double w = (S * c * p + 1.0) / (c * Tr);
    return (-Q / (2.0 * M_PI * Tr * p) * std::cyl_bessel_k(0.0, r * sqrt(w)));
  };
};

struct BarkerHerbert
{
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
                t_1(t_1), t_2(t_2), s_1(s_1), s_2(s_2), Q(Q), prec(prec)
  {
    tau = std::log(2.0) / time.array();
  };
  double lp(double p)
  {

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

    double term_1 = (bk_n1 * bk_a0 - bk_a1 * bk_n0 * ct);
    double term_2 = (bk_n0 * bi_n1 + bk_n1 * bi_n0);

    double drawdown;

    // if radius is inside the patch
    if (radius <= radius_patch)
    {

      drawdown = std::cyl_bessel_k(0.0, N * radius) / p +
        (term_1 * std::cyl_bessel_i(0.0, N * radius)) / denom;

      return (drawdown);
    }

    // if radius is outside the patch
    drawdown = term_2 * std::cyl_bessel_k(0.0, A * radius) / denom;

    return (drawdown * Q / (2.0 * M_PI * t_1));
  }
};

// [[Rcpp::export]]
Eigen::RowVectorXd stehfest_v(
    int n)
{

  assert(n < 20);
  assert(n % 2 == 0);

  int n_div_2 = n / 2;
  int s, e;
  double z;

  Eigen::RowVectorXd fact(n + 1);
  Eigen::RowVectorXd v(n);

  // calculate factorials
  for (unsigned int i = 0; i <= n; ++i)
  {
    fact[i] = std::tgamma(i + 1);
  }

  for (int i = 0; i < n; ++i)
  {
    if (i > n_div_2 - 1)
    {
      s = n_div_2 - 1;
    }
    else
    {
      s = i;
    }
    e = (int)(i / 2);

    z = 0.0;
    for (int k = s; k >= e; --k)
    {
      z += (pow((double)k + 1.0, n_div_2) * fact[2 * k + 2]) /
        (fact[n_div_2 - k - 1] * fact[k + 1] * fact[k] *
          fact[i - k] * fact[2 * k - i + 1]);
    }

    v[i] = pow(-1.0, n_div_2 + i + 1.0) * z;
  }

  return (v);
}

// [[Rcpp::export]]
Eigen::MatrixXd stehfest_p(
    Eigen::VectorXd time,
    int n_terms)
{

  Eigen::RowVectorXd v = Eigen::RowVectorXd::LinSpaced(n_terms, 1.0, (double)n_terms);

  return (time * v);
}


// // [[Rcpp::export]]
// std::complex<double> tth(
//     std::complex<double> time)
// {
//   return(std::tanh(time));
// }
//
// // [[Rcpp::export]]
// std::complex<double> cch(
//     std::complex<double> time)
// {
//   return(std::cosh(time));
// }

// // [[Rcpp::export]]
// Rcpp::ComplexVector bbl(
//     Rcpp::ComplexVector time,
//     double nu,
//     bool expon_scaled,
//     size_t n_seq)
// {
//   return(specialfunctions::bessel_k_cplx(time, 0, false, 1));
// }


template <typename T>
Eigen::VectorXd stehfest(T &well, int n_terms)
{

  Eigen::VectorXd v = stehfest_v(n_terms);
  Eigen::MatrixXd p = stehfest_p(well.tau, n_terms);

  std::vector<double> p_vec(p.data(), p.data() + p.size());

  for (auto &out : p_vec)
    out = well.lp(out);

  p = Eigen::Map<Eigen::MatrixXd>(p_vec.data(), p.rows(), p.cols());
  Eigen::VectorXd ret = (p * v).array() * (well.tau).array();

  return(ret);
}

template <typename T>
Eigen::VectorXd stehfest_xy(T &well, int n_terms)
{

  Eigen::VectorXd v = stehfest_v(n_terms);
  Eigen::MatrixXd p = stehfest_p(well.tau, n_terms);

  // std::vector<double> p_vec(p.data(), p.data() + p.size());
  unsigned int k = 0;
  double z = 0.0;
  double x = 0.0;
  double t_p = 0.0;

  for (unsigned int i = 0; i < p.rows(); ++i)
  {
    z = well.z[i];
    x = well.x[i];
    t_p = well.time[i];
    for (unsigned int j = 0; j < p.cols(); ++j)
    {
      p(i, j) = well.lp_xy(p(i, j), t_p, z, x);
      k++;
    }
  }

  // p = Eigen::Map<Eigen::MatrixXd>(p_vec.data(), p.rows(), p.cols());
  Eigen::VectorXd ret = (p * v).array() * (well.tau).array();

  return(ret);
}

// https://mpmath.org/doc/current/calculus/inverselaplace.html#mpmath.calculus.inverselaplace.FixedTalbot.calc_laplace_parameter
// [[Rcpp::export]]
Eigen::MatrixXcd cohen_p(
    Eigen::VectorXd time,
    unsigned int n_terms)
{

  unsigned int n_t = time.size();

  Eigen::VectorXd time_2 = 2.0 * time;

  std::complex<double> j(0.0, 1.0);


  Eigen::VectorXd gamma = (2.0 / 3.0) *
    (std::floor(n_terms / 1.31) + log(10.0) + (time_2.array().log()));

  // 0 to n_terms
  Eigen::RowVectorXcd k = Eigen::RowVectorXcd::LinSpaced(n_terms + 1, 0.0, (double)n_terms) * j * M_PI;

  Eigen::VectorXcd term_1 = Eigen::VectorXcd::Zero(n_t);

  term_1.real() = gamma.array() / (time_2.array());

  Eigen::MatrixXcd term_2 = (time.array().inverse()).matrix() * k;

  term_2.colwise() += term_1;

  return (term_2);
}


// [[Rcpp::export]]
Eigen::VectorXd cohen_c(
    double d,
    unsigned int n_terms)
{
  double n = double(n_terms);
  // double m = n_terms + 1.0;
  double c_dbl = -d;
  double b = -1.0;
  Eigen::VectorXd c(n_terms);

  for (unsigned int k = 0; k < n_terms; ++k) {
    c_dbl = (double)b - c_dbl;
    b = 2.0 * ((double)k + n) * ((double)k - n) * b /
      ((2.0 * (double)k + 1.0) * ((double)k + 1.0));

    c(k) = c_dbl;
  }

  return (c);
}


template <typename T>
Eigen::VectorXd cohen_xy(T &well, unsigned int n_terms)
{

  Eigen::MatrixXcd p = cohen_p(well.time, n_terms);
  Eigen::VectorXd ret(well.time.size());

  double z = 0.0;   // distance along fracture
  double x = 0.0;   // distance into matrix
  double t_p = 0.0; // time at p
  double f_p = 0.0; // laplace solution for p


  // acceleration variables
  double s;
  double d = pow(3.0 + sqrt(8), n_terms);
  d = (d + 1.0 / d) / 2.0;
  Eigen::VectorXd c = cohen_c(d, n_terms);


  for (unsigned int i = 0; i < p.rows(); ++i)
  {
    z = well.z[i];
    x = well.x[i];
    t_p = well.time[i];

    // this is the first term A[0]
    ret[i] = well.cohen_xy(p(i, 0), t_p, z, x);


    // this is the summation term
    // need to figure out how to large values in sinh so that oscillations
    // do not occur
    s = 0.0;
    for (unsigned int j = 1; j < (n_terms + 1); ++j)
    {
      f_p = well.cohen_xy(p(i, j), t_p, z, x);
      s += c[j-1] * f_p;
    }

    ret[i] = (std::exp(well.gamma[i] / 2.0) / t_p) * ((ret[i] / 2.0) - (s / d));
  }


  // ret = ((well.gamma / 2.0).array().exp() / well.time.array()) * ret.array();

  return(ret);
}

// template <typename T>
// Eigen::VectorXd cohen(T &well, unsigned int n_terms)
// {
//
//   unsigned int n_times = well.time.size();
//   Eigen::MatrixXcd p = cohen_p(well.time, n_terms);
//
//   // acceleration variables
//   double d = pow(3.0 + sqrt(8), n_terms);
//   d = (d + 1.0 / d) / 2.0;
//   Eigen::VectorXd c = cohen_c(d, n_terms);
//
//   for (auto &f_p : p)
//     f_p = well.cohen(f_p);
//
//   Eigen::VectorXd ret = p.col(0);
//   Eigen::VectorXd s = p.rightCols(n_terms) * c;
//
//   ret = (well.gamma.array() / 2.0).exp() / well.time.array() *
//     ((ret.array() / 2.0) - (s.array() / d));
//
//   return(ret);
// }


// [[Rcpp::export]]
Eigen::VectorXd cooper_bredehoeft_papadopulos_laplace(
    Eigen::VectorXd time,
    double r,
    double r_c,
    double r_w,
    double Tr,
    double S,
    double h_0,
    int n_terms)
{

  Rcpp::Rcout << "The value is r " << r << std::endl;
  Rcpp::Rcout << "The value is r_c " << r_c << std::endl;
  Rcpp::Rcout << "The value is r_w " << r_w << std::endl;
  Rcpp::Rcout << "The value is Tr " << Tr << std::endl;
  Rcpp::Rcout << "The value is S " << S << std::endl;
  Rcpp::Rcout << "The value is h_0 " << h_0 << std::endl;


  CooperBredehoeftPapadopulos well(time, r, r_c, r_w, Tr, S, h_0);
  Eigen::VectorXd out = stehfest(well, n_terms);

  // time equal to zero replace with h_0
  out = out.unaryExpr([h_0](double v) { return std::isfinite(v)? v : h_0; });
  return (out);// / (2.0 * M_PI * Tr));
}


// [[Rcpp::export]]
Eigen::VectorXd papadopulos_cooper_laplace(
    Eigen::VectorXd time,
    double Q,
    double r,
    double r_c,
    double r_w,
    double Tr,
    double S,
    double prec,
    int n_terms)
{

  PapadopulosCooper well(time, Q, r, r_c, r_w, Tr, S, prec);

  return (stehfest(well, n_terms));// / (2.0 * M_PI * Tr));
}

// [[Rcpp::export]]
Eigen::VectorXd jacob_lohman_laplace(
    Eigen::VectorXd time,
    double s,
    double r,
    double Tr,
    double S,
    double prec,
    int n_terms)
{

  JacobLohman well(time, s, r, Tr, S, prec);

  return (stehfest(well, n_terms) * (2.0 * M_PI * r * Tr));
}

// [[Rcpp::export]]
Eigen::VectorXd hantush_jacob_laplace(
    Eigen::VectorXd time,
    double c,
    double r,
    double Tr,
    double S,
    double Q,
    double prec,
    int n_terms)
{

  HantushJacob well(time, c, r, Tr, S, Q, prec);

  return (stehfest(well, n_terms));
}

// [[Rcpp::export]]
Eigen::VectorXd barker_herbert(
    Eigen::VectorXd time,
    double radius,
    double radius_patch,
    double t_1,
    double t_2,
    double s_1,
    double s_2,
    double Q,
    double prec,
    int n_terms)
{

  BarkerHerbert well(time, radius, radius_patch,
                     t_1, t_2, s_1, s_2, Q, prec);

  return (stehfest(well, n_terms));
}


// [[Rcpp::export]]
Eigen::VectorXd parallel_fractures_solute(
    Eigen::VectorXd time,
    Eigen::VectorXd z,
    Eigen::VectorXd x,
    Eigen::VectorXd concentration_influent,
    Eigen::VectorXd time_influent,
    double c_0,
    double b,
    double B,
    double v,
    double alpha_l,
    double D_star,
    double k_f,
    double k_m,
    double t_half,
    double rho_b,
    double theta,
    double tortuosity, // τ
    unsigned int n_terms
)
{
  ParallelFracturesSolute fracture(
      time,
      z,
      x,
      concentration_influent,
      time_influent,
      c_0,
      b,
      B,
      v,
      alpha_l,
      D_star,
      k_f,
      k_m,
      t_half,
      rho_b,
      theta,
      tortuosity, // τ
      n_terms
  );

  return (cohen_xy(fracture, n_terms));
}


// [[Rcpp::export]]
Eigen::VectorXd parallel_fractures_heat(
    Eigen::VectorXd time,
    Eigen::VectorXd z,
    Eigen::VectorXd x,
    Eigen::VectorXd temperature_influent,
    Eigen::VectorXd time_influent,
    double t_0,
    double b,
    double B,
    double v,
    double lambda_fracture, // thermal conductivity of water
    double lambda_matrix,   // thermal conductivity of solids
    double spec_heat_w,     // specific heat of water
    double spec_heat_s,     // specific heat of solids
    double rho_w,           // density of water
    double rho_s,           // density of solids
    double theta,
    unsigned int n_terms)
{


  ParallelFracturesHeat heat(time,
                             z,
                             x,
                             temperature_influent,
                             time_influent,
                             t_0,
                             b,
                             B,
                             v,
                             lambda_fracture, // thermal conductivity of water
                             lambda_matrix,   // thermal conductivity of solids
                             spec_heat_w,     // specific heat of water
                             spec_heat_s,     // specific heat of solids
                             rho_w,           // density of water
                             rho_s,           // density of solids
                             theta,
                             n_terms);

  return (cohen_xy(heat, n_terms));
}




// // [[Rcpp::export]]
// Rcpp::List lap2(Eigen::MatrixXd p,
//                 Eigen::VectorXd tau,
//                 Eigen::VectorXd v,
//                 Eigen::VectorXd z,
//                 Eigen::VectorXd x,
//                 double lambda,
//                 double nu,
//                 double sigma,
//                 double c_0,
//                 double a,
//                 double k2,
//                 double B,
//                 double g)
// {

//   unsigned int n_p_r = p.rows();
//   unsigned int n_p_c = p.cols();
//   unsigned int n_z = z.size();
//   unsigned int n_x = x.size();

//   double p, sqrt_p;

//   Eigen::VectorXd term_z_pos = (nu * z).exp();    // fracture terms
//   Eigen::VectorXd term_z_neg = (-nu * z).exp();   // fracture terms
//   Eigen::VectorXd term_x = (0.5 * B) - x.array(); // matrix term

//   Eigen::VectorXd fracture(n_z);
//   Eigen::VectorXd matrix(n_x);

//   Eigen::MatrixXd to_sum(n_z, n_x);
//   Rcpp::List ret(n_p_r);

//   for (unsigned int i = 0; i < n_p_r; ++i)

//     for (unsigned int j = 0; j < n_p_c; ++i)
//     {
//       to_sum.setZero();
//       {
//         sqrt_p = sqrt(p[i, j]);

//         fracture = term_z_pos * (c_0 / (p[i] - lambda));
//         fracture.array() /= (term_z_neg.array() *
//                              (sqrt(1.0 + k2 * (p[i] + sqrt(p[i] / a * std::tanh(sigma * sqrt_p))))));

//         matrix = (term_x * (g * sqrt_p)).cosh() / std::cosh(sigma * sqrt_p);
//         to_sum += (fracture * matrix.transpose()) * v[i];
//       }

//       ret[i] = to_sum * tau[i];
//     }

//   return (ret);
// }

/*** R
n <- 10000
time = c(0, 1:86400)
# CooperBredehoeftPapadopulos well(time, r, r_c, r_w, Tr, S, h_0);

kern_slug <- frecipes:::cooper_bredehoeft_papadopulos_laplace(time,
                        r = 0.10,
                        r_c = 0.10,
                        r_w = 0.10,
                        S = 1e-5,
                        Tr = 5e-4,
                        h_0 = 1,
                        n = 14L)


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
  a <- frecipes:::jacob_lohman_laplace(times, rw, Tr, s, S, prec, n_terms),
  b <- frecipes:::hantush_jacob(times, flow_rate, rw,S,Tr, lab, prec),
  # c <- frecipes:::hantush_jacob_quad(times, lab, rw, Tr, S, Q, 1e-16),
  # c <- frecipes:::barker_herbert(times, c, rw, Tr, S, Q, prec, 12L),
  check = FALSE
)

bench::mark(
  a <- frecipes:::hantush_jacob_laplace(times, cc, rw, Tr, S, Q, prec, n_terms),
  b <- frecipes:::hantush_jacob(times, flow_rate, rw,S,Tr, lab, prec),
  # c <- frecipes:::hantush_jacob_quad(times, lab, rw, Tr, S, Q, 1e-16),
  # c <- frecipes:::barker_herbert(times, c, rw, Tr, S, Q, prec, 12L),
  check = FALSE
)



times <- 10^seq(-5, 2, by = 0.1)
s <- 10
Tr = 10 # transmissivity of aquifer, m^2/d
S = 1e-5 # storage coefficient of aquifer, -
rw = 0.15 # radius of well, m
a <- frecipes:::jacob_lohman_laplace(times, rw, Tr, s, S, prec, 16L)[[1]]
a <- (2.0 * pi * rw * Tr) * a
plot(y = abs(a), x=times, type = "l", log = "xy")



times <- 10^seq(-5, 2, by = 0.1)
Tr = 10 # transmissivity of aquifer, m^2/d
S = 1e-4 # storage coefficient of aquifer, -
rw = 0.15 # radius of well, m
rc = 0.15
r = 0.15
Q = 10
a <- frecipes:::papadopulos_cooper_laplace(times,
                                           Q,
                                           r,
                                           rc,
                                           rw,
                                           Tr,
                                           S,
                                           prec, 16L)

Eigen::VectorXd time,
double Q,
double r,
double r_c,
double r_w,
double Tr,
double S,
double prec,
int n_terms
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
