#include "hydrorecipes.h"

// from WTAQ manual pg 11
// [[Rcpp::export]]
double well_bore_storage(double rc,  // casing radius where water level is changing
                         double rw,  // well screen radius
                         double Ss,  // specific storage
                         double zpl, // is the depth below top of aquifer or initial water table to the bottom of the screened interval of the pumped well (units of length);
                         double zpd) // is the depth below top of aquifer or initial water table to the top of the screened interval of the pumped well (units of length). [Note that is the length of the screened interval of the pumped well.]

{
  return(M_PI * rc * rc / (2.0 * M_PI * rw * rw * Ss * (zpl - zpd)));
}

// from WTAQ manual pg 11
// [[Rcpp::export]]
double well_skin(double kr, // aquifer radial hydraulic conductivity
                 double ks, // is the hydraulic conductivity of the well-bore skin (units of length per time); and
                 double ds, // is the thickness of the well-bore skin (units of length);
                 double rw) // well screen radius
{
  return((kr * ds) / (ks * rw));
}


// from WTAQ manual pg 12
// [[Rcpp::export]]
double shape_factor(double L,
                    double rp,
                    double kr,
                    double kz)
{
  double m = sqrt(kr / kz);
  double x = m * L / (2.0 * rp);
  return(L / (log(x + std::sqrt(1 + x * x))));
}

// from WTAQ manual pg 12
// [[Rcpp::export]]
double well_delay(double rp,
                  double rw,
                  double Ss,
                  double kr,
                  double kz,
                  double L)
{
  double fp = shape_factor(L, rp, kr, kz);
  return( M_PI * rp * rp / (2.0 * M_PI * rw * rw * Ss * fp));
}




// Moench, A.F., 1997. Flow to a well of finite diameter in a homogeneous, anisotropic water table aquifer. Water Resources Research, 33(6), pp.1397-1407.
// Table 1
// [[Rcpp::export]]
double dimensionless_time(double Kr,
                          double t,
                          double rw,
                          double Ss,
                          double b = 1.0) {
  return ((Kr * t * b) / (rw * rw * Ss));
}

// [[Rcpp::export]]
double dimensionless_head(double Kr,
                          double t,
                          double rw,
                          double ht,
                          double h,
                          double Q,
                          double b = 1) {
  return ((4.0 * M_PI * Kr * b * (ht - h) / Q));
}


// r/rw
// rw/b
// z/b
// Kz/Kr
// l/b
// d/b
// [[Rcpp::export]]
double dimensionless_ratio(double a,
                           double b) {
  return (a / b);
}


// Kd, rwd
// Beta, rd
// [[Rcpp::export]]
double dimensionless_beta(double Kd,
                          double rd) {
  return (Kd * rd * rd);
}


// [[Rcpp::export]]
double dimensionless_sigma(double Ss,
                           double b,
                           double Sy) {
  return (Ss * b / Sy);
}

// [[Rcpp::export]]
double dimensionless_alpha(double Ss,
                           double b,
                           double Sy) {
  return (Ss * b / Sy);
}

// [[Rcpp::export]]
double dimensionless_gamma(double alpha,
                           double b,
                           double Sy,
                           double Kz) {
  return (alpha * b * Sy / Kz);
}
// [[Rcpp::export]]
double dimensionless_well_bore_storage(double rc,
                       double rw,
                       double Ss,
                       double l,
                       double d) {
  return (M_PI * rc * rc / (2.0 * M_PI * rw * rw * Ss * (l - d)));
}
// [[Rcpp::export]]
double dimensionless_w_prime(double rc,
                             double rw,
                             double Ss,
                             double f_prime) {
  return (M_PI* rc*rc / (2.0 * M_PI * rw * rw * Ss * f_prime));
}

// [[Rcpp::export]]
double dimensionless_S(double Kr,
                       double ds,
                       double Ks,
                       double rw) {
  return (Kr * ds / (Ks * rw));
}


/*** R




*/
