#include "hydrorecipes.h"

//==============================================================================
//' @title
//' well_function_coefficient
//'
//' @description
//' Calculation Coefficient Q/(4 pi * T)
//'
//' @param flow_rate well flow rates
//' @param transmissivity aquifer transmissivity
//'
//' @return coefficient for Theis and Hantush well function
//'
//'
//' @export
//'
// [[Rcpp::export]]
double well_function_coefficient(const double flow_rate,
                                 const double transmissivity)
{
  return(flow_rate / (4.0 * M_PI * transmissivity));
}

//==============================================================================
//' @title
//' grf_coefficient
//'
//' @description
//' Coefficient for the grf without pumping
//'
//' @param flow_rate (vector) the flow rate
//' @param radius (double) distance to center of source
//' @param K (double) hydraulic conductivity of the fracture system
//' @param thickness (double) thickness
//' @param flow_dimension (double) flow dimension
//'
//' @return coefficient
//'
//'
//' @export
//'
// [[Rcpp::export]]
double grf_coefficient(const double flow_rate,
                       const double radius,
                       const double K,
                       const double thickness,
                       const double flow_dimension)
{

  double v = 1.0 - flow_dimension / 2.0;

  return ((flow_rate * pow(radius, (2.0 * v))) /
          (4.0 * pow(M_PI, (1.0 - v)) * K *
           pow(thickness, (3.0 - flow_dimension))));
}

//==============================================================================
//' @title
//' hantush_epsilon
//'
//' @description
//' Calculation of r^2/(4B^2)
//'
//' @param radius distance to monitoring well
//' @param leakage aquifer transmissivity
//'
//' @return coefficient Hantush well function
//'
//'
//' @export
//'
// [[Rcpp::export]]
double hantush_epsilon(const double radius,
                       const double leakage)
{

  return((radius * radius) / (4.0 * leakage * leakage));
}

//==============================================================================
//' @title
//' grf_u_time
//'
//' @description
//' Calculation of grf u
//'
//' @param radius distance to monitoring interval
//' @param storativity aquifer storativity
//' @param K aquifer hydraulic conductivity
//' @param time prediction times
//'
//' @return u for well function
//'
//'
//' @export
//'
// [[Rcpp::export]]
double grf_u(const double radius,
             const double storativity,
             const double K,
             const double time)
{

  return ((radius * radius * storativity) /
         (4.0 * K * time));
  ;
}
