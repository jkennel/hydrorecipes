#include "frecipes.h"

// [[Rcpp::export]]
double any_decimal(std::vector<double> x)
{

  double div = 1.0;
  double prev = 1.0;

  for (unsigned int i = 0; i < x.size(); ++i) {
    div = std::fmod(abs(x[i]), 1.0);

    if(div > 0.0) {
      prev = std::min(div, prev);
    }

  }

  if(prev == 1.0) {
    return(1.0);
  }
  if(prev >= 0.1) {
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
    if(out == 1) {
      return(1);
    }
  }
  return(out);
}




/*** R

bench::mark(
  frecipes:::gcd(sample(seq(0, 1000000, 10), 20, replace = TRUE)),
  frecipes:::any_decimal(c(sample(seq(0, 1000000, 10)), 1.1, 2.02)),
  check = FALSE
)

*/
