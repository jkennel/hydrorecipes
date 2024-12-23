#include "hydrorecipes.h"

// [[Rcpp::export]]
Eigen::VectorXd row_sums_eigen(Eigen::Map<Eigen::MatrixXd>& x) {

  if (x.cols() == 1) {
    return(x);
  }

  return(x.rowwise().sum());

}

// [[Rcpp::export]]
Rcpp::NumericVector row_multiply3(Rcpp::List x,
                                 Rcpp::NumericVector y) {

  Rcpp::NumericVector z = x[0];
  Rcpp::NumericVector m = z;
  z = z * y[0];

  for (int i = 1; i < x.size(); ++i) {
    m = x[i];
    z = z + m * y[i];
  }

  return(z);
}

// [[Rcpp::export]]
Eigen::VectorXd row_multiply2(Rcpp::List x,
                              Eigen::VectorXd y) {

  Eigen::VectorXd x_col = Rcpp::as<Eigen::Map<Eigen::VectorXd>>(x[0]) * y[0];

  for (int i = 1; i < x.size(); ++i) {
    Eigen::VectorXd m = Rcpp::as<Eigen::Map<Eigen::VectorXd>>(x[i]) * y[i];
    x_col += m;
  }

  return(x_col);
}

// [[Rcpp::export]]
Eigen::VectorXd list_multiply_subset(const Rcpp::List x,
                              const Eigen::VectorXd& y,
                              const Eigen::VectorXi& ind) {

  Eigen::VectorXd x_col = Rcpp::as<Eigen::Map<Eigen::VectorXd>>(x[0])(ind) * y[0];

  for (int i = 1; i < x.size(); ++i) {
    Eigen::VectorXd m = Rcpp::as<Eigen::Map<Eigen::VectorXd>>(x[i])(ind) * y[i];
    x_col += m;
  }

  return(x_col);
}



// [[Rcpp::export]]
double any_decimal(std::vector<double> x)
{

  double div = 1.0;
  double prev = 1.0;

  for (unsigned int i = 0; i < x.size(); ++i) {
    div = std::fmod(abs(x[i]), 1.0);

    if (div > 0.0) {
      prev = std::min(div, prev);
    }

  }

  if (prev == 1.0) {
    return(1.0);
  }
  if (prev >= 0.1) {
    return(10.0);
  }
  if (prev >= 0.01) {
    return(100.0);
  }
  if (prev >= 0.001) {
    return(1000.0);
  }
  if (prev >= 0.0001) {
    return(10000.0);
  }

  return (1.0);
}

// [[Rcpp::export]]
std::vector<double> decimal_to_scaled_integer(std::vector<double> x)
{
  double mult = any_decimal(x);

  if (mult != 1.0) {
    for (auto &out : x)
      out = out * mult;
  }

  return(x);
}

// [[Rcpp::export]]
unsigned int gcd(std::vector<unsigned int> x)
{
  unsigned int out = x[0];
  for (unsigned int i = 0; i < x.size(); ++i) {
    out = std::gcd(out, x[i]);
    if (out == 1) {
      return(1);
    }
  }
  return(out);
}




/*** R

bench::mark(
  hydrorecipes:::gcd(sample(seq(0, 1000000, 10), 20, replace = TRUE)),
  hydrorecipes:::any_decimal(c(sample(seq(0, 1000000, 10)), 1.1, 2.02)),
  check = FALSE
)

*/
