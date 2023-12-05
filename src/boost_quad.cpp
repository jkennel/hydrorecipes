#include "frecipes.h"

using namespace boost::math::quadrature;

// [[Rcpp::export]]
double hantush_jacob_gauss_kronrod(double t,
            double lab,
            double r,
            double T,
            double S,
            double Q,
            double prec) {


  // eq 14 Hantush and Jacob 1955 (short times)
  // eq 8 Hantush and Jacob 1955
  auto f1 = [&](double tau) {
    return (1.0 / tau) * std::exp(-tau - r * r / (4.0 * lab * lab * tau));};
  double u = (S * r * r) / (4.0 * T * t);


  double error;
  double out = gauss_kronrod<double, 15>::integrate(f1,
                                                    u,
                                                    std::numeric_limits<double>::infinity(),
                                                    5,
                                                    prec);

  return(-Q / (4.0 * M_PI * T) * out);
}

// [[Rcpp::export]]
std::vector<double> hantush_jacob_quad(std::vector<double> t,
            double lab,
            double r,
            double T,
            double S,
            double Q,
            double prec) {

    for (auto &out : t)
    out = hantush_jacob_gauss_kronrod(out, lab, r, T, S, Q, prec);

  return(t);
}



/*** R

*/
