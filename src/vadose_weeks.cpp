#include "frecipes.h"

//==============================================================================
// [[Rcpp::export]]
double weeks_1979(const double lag,
                  const double D,
                  const double L,
                  const double precision,
                  const bool inverse) {

  const double d_term = M_2PI * D * lag / (4.0 * L * L);
  double term_val = 0.0;
  double ret = 0.0;
  double exp_val = 0.0;
  unsigned int m = 1;
  bool more_precise = TRUE;

  if (inverse) {
    if (d_term < 0.001) {
      return(0.0);
    }
  } else {
    if (d_term < 0.001) {
      return(1.0);
    }
  }


  while(more_precise) {
    exp_val  = -double(m * m) * d_term;
    term_val = std::pow(-1.0, (double(m) - 1.0) / 2.0) / double(m) * std::exp(exp_val);
    ret += term_val;
    m   += 2;
    more_precise = std::abs(term_val) > precision;
  }

  if (inverse) {
    ret = (1.0 - (4.0 / M_PI) * ret);
    if (ret < 0.0) {
      ret = 0.0;
    }
  } else {
    ret = ((4.0 / M_PI) * ret);
    if (ret > 1.0) {
      ret = 1.0;
    }
  }


  return(ret);

}

//==============================================================================
//' @title
//' vadose_response
//'
//' @description
//' weeks_1979 1-D air diffusivity
//'
//' @param time numeric vector of elapsed times
//' @param air_diffusivity A numeric value of the unsaturated zone air diffusivity
//' @param thickness A numeric value of the unsaturated zone thickness
//' @param precision A numeric value of for the solution precision
//' @param inverse A logical value indicating if an inverse water level relationship is desired
//'
//' @return weeks 1979 model
//'
//'
//' @export
//'
//' @examples
//' vr <- vadose_response(time = 0:43200,
//'                        air_diffusivity = 0.20,
//'                        thickness = 40,
//'                        precision = 1e-10,
//'                        inverse = FALSE)
//'
//' @noRd
//'
// [[Rcpp::export]]
Rcpp::List vadose_response(std::vector<double> time,
                                     const double air_diffusivity,
                                     const double thickness,
                                     const double precision,
                                     const bool inverse) {

  for (auto &out : time)
    out =  weeks_1979(
      out,
      air_diffusivity,
      thickness,
      precision,
      inverse);

  return Rcpp::List::create(Rcpp::Named("vadose_weeks") = time);


  // return(time);
}


// [[Rcpp::export]]
Rcpp::NumericVector vadose_response2(const Rcpp::NumericVector time,
                                     double air_diffusivity,
                                     double thickness,
                                     double precision,
                                     bool inverse) {

  unsigned int n = time.size();

  Rcpp::NumericVector output(n);

  for (unsigned int i = 0; i < n; ++i) {
    output(i) = weeks_1979(
      time(i),
      air_diffusivity,
      thickness,
      precision,
      inverse);
  }


  return(output);
}





/***R

bench::mark(
tmp <- frecipes:::vadose_response(time = as.numeric(0:(43200*100)),
                       air_diffusivity = 0.20, thickness = 40,
                       precision = 1e-12,
                       inverse = FALSE),
# tmp2 <- frecipes:::vadose_response2(time = as.numeric(0:(43200*100)),
#                        air_diffusivity = 0.20, thickness = 40,
#                        precision = 1e-12,
#                        inverse = FALSE),
relative = FALSE
)


*/
