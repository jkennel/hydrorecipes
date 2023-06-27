// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppThread)]]

#define EIGEN_FFTW_DEFAULT

#include <RcppEigen.h>
#include <unsupported/Eigen/FFT>
#include <unsupported/Eigen/Splines>
#include <fftw3.h>
#include <Eigen/StdVector>
#include <specialfunctions.h>
#include <RcppThread.h>

using namespace Eigen;
using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::MatrixXcd;
using Eigen::Vector2d;
using Eigen::Vector3d;
using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::VectorXcd;
using Eigen::RowVectorXd;
using Eigen::ArrayXd;
using Eigen::ArrayXi;

using Eigen::FFT;

using namespace Rcpp;


//==============================================================================
// to_dummy.cpp
IntegerMatrix to_dummy(const IntegerVector& x,
                       size_t n_fact);


//==============================================================================
// lags.cpp
int check_lag(int n, int lag, int n_shift);
int get_length(int n, int n_subset);
int get_start(int n_out, int lag, int n_subset);
int get_end(int n, int n_out, int lag, int n_subset);
Eigen::ArrayXd log_lags_eigen(size_t n, size_t max_lag);
Eigen::VectorXi log_lags_eigen2(size_t n, size_t max_lag);

NumericVector shift_subset(const NumericVector& x,
                           int lag = 0,
                           int n_subset = 1,
                           int n_shift = 0);
NumericMatrix lag_matrix(const NumericMatrix& x,
                         const IntegerVector& lags,
                         CharacterVector suffix,
                         std::string prefix,
                         int n_subset = 1,
                         int n_shift = 0);
Eigen::MatrixXd distributed_lag_parallel(const Eigen::VectorXd& x,
                                         const Eigen::MatrixXd& bl,
                                         int lag_max,
                                         int n_subset = 1,
                                         int n_shift = 0);

//==============================================================================
// harmonic.cpp
Eigen::MatrixXd harmonic_double(const Eigen::VectorXd& time,
                                const Eigen::RowVectorXd& frequency,
                                double cycle_size);

//==============================================================================
// fft_helpers.cpp
size_t index_from_i_j(size_t i, size_t j, size_t n_col);
size_t index_from_j_i(size_t i, size_t j, size_t n_col);

size_t get_column_number(size_t n);

size_t next_n_eigen(size_t n);
Eigen::VectorXd pad_vector(Eigen::VectorXd x, size_t n_old, size_t n_new);

Eigen::MatrixXd detrend_matrix(const Eigen::MatrixXd& x);
Eigen::MatrixXd demean_matrix(const Eigen::MatrixXd& x);
Eigen::MatrixXd detrend_and_demean_matrix(const Eigen::MatrixXd& x,
                                          bool detrend,
                                          bool demean);

Eigen::VectorXd modified_daniell(Eigen::VectorXi spans);
Eigen::MatrixXcd kernel_apply(Eigen::MatrixXcd& x,
                              Eigen::VectorXd& y);
Eigen::VectorXd spec_taper(size_t n_row, double p = 0.1);



Eigen::VectorXi make_groups(size_t n_groups,
                            size_t n);
Eigen::ArrayXd power_spaced(size_t n, double min, double max, double power);
Eigen::VectorXd group_frequency(Eigen::ArrayXd frequencies,
                                size_t n_groups);
Eigen::VectorXd determine_frequency(size_t n);
Eigen::MatrixXcd check_ffts(Eigen::MatrixXcd& x,
                            double cutoff);
Eigen::VectorXi which_indices(const Eigen::VectorXd& x,
                              const Eigen::VectorXd& knots);



//==============================================================================
// b_spline.cpp
Eigen::MatrixXd b_spline(const Eigen::ArrayXd& x,
                         const Eigen::ArrayXd& knots,
                         size_t degree = 3);


//==============================================================================
//fft_windows.cpp
Eigen::VectorXd window_hann(size_t n);
Eigen::VectorXcd window_hann_cplx(size_t n);
Eigen::VectorXd window_rectangle(size_t n);
Eigen::VectorXd window_tukey(size_t n, double r);
double window_scale(Eigen::VectorXd window, size_t n_new, size_t n_fft);



//==============================================================================
//fft.cpp
Eigen::MatrixXcd fft_matrix(Eigen::MatrixXd x,
                            size_t n_new);

Eigen::MatrixXcd multiply_ffts(Eigen::MatrixXcd& x);

//convolve
Eigen::VectorXd convolve_vec(Eigen::VectorXd x, Eigen::VectorXd y);
Eigen::VectorXd convolve_overlap_add(Eigen::VectorXd x, Eigen::VectorXd y);
Eigen::VectorXd convolve_overlap_save(Eigen::VectorXd x, Eigen::VectorXd y);
Eigen::VectorXd convolve_tf(Eigen::VectorXd x, Eigen::VectorXcd y);
Eigen::MatrixXd convolve_matrix(const Eigen::VectorXd& x,
                                const Eigen::MatrixXd& y,
                                bool remove_partial = true,
                                bool reverse = true);
Eigen::VectorXd convolve_filter(const Eigen::VectorXd& x,
                               const Eigen::VectorXd& y,
                               bool remove_partial,
                               bool reverse);


// Spectrum
Eigen::MatrixXcd spec_welch(Eigen::MatrixXd& x,
                            size_t length_subset,
                            double overlap,
                            Eigen::VectorXd window = Eigen::VectorXd::Zero(0));
// Eigen::MatrixXcd spec_welch_trunc(const Eigen::MatrixXd& x,
//                                   size_t length_subset,
//                                   double overlap,
//                                   Eigen::VectorXd window = Eigen::VectorXd::Zero(0));
// Eigen::MatrixXcd spec_welch_complete(const Eigen::MatrixXd& x,
//                                      size_t length_subset,
//                                      double overlap,
//                                      Eigen::VectorXd window = Eigen::VectorXd::Zero(0));


Eigen::MatrixXcd spec_pgram(Eigen::MatrixXd& x,
                            const Eigen::VectorXi& spans,
                            bool detrend,
                            bool demean,
                            double taper = 0.1);


// Solve

Eigen::MatrixXcd solve_cplx_parallel(const Eigen::MatrixXcd& x);
Eigen::MatrixXcd solve_cplx_irr(Eigen::MatrixXcd& x,
                                size_t n_groups);

Eigen::MatrixXcd transfer_pgram_smooth(Eigen::MatrixXd& x,
                                       const Eigen::VectorXi& spans,
                                       bool detrend,
                                       bool demean,
                                       double taper,
                                       double power,
                                       size_t n_groups);
Eigen::MatrixXcd transfer_pgram(Eigen::MatrixXd& x,
                                const Eigen::VectorXi& spans,
                                bool detrend,
                                bool demean,
                                double taper,
                                double power,
                                size_t n_groups);

Eigen::MatrixXcd transfer_welch(Eigen::MatrixXd& x,
                                size_t length_subset,
                                double overlap,
                                Eigen::VectorXd window = Eigen::VectorXd::Zero(0));



// Processing
Eigen::MatrixXd ordinary_coherence_phase(const Eigen::MatrixXcd& x);
Eigen::MatrixXd frequency_to_time_domain(Eigen::MatrixXcd& pgram,
                                         size_t n_groups);
// Eigen::ArrayXd frf_to_brf(const Eigen::VectorXcd& x,
//                           std::complex<double> dc1,
//                           std::complex<double> dc2);
Eigen::MatrixXcd interpolate_tf(Eigen::MatrixXcd& x,
                                 const Eigen::ArrayXd& frequency_irregular,
                                 const Eigen::ArrayXd& frequency_regular,
                                 Eigen::VectorXd& knots);
// Eigen::MatrixXd predict_pgram_frf(Eigen::MatrixXd& x,
//                                   Eigen::VectorXi span);
Eigen::MatrixXd predict_pgram_frf(Eigen::MatrixXd& x,
                                  Eigen::MatrixXd& x_out,
                                  Eigen::VectorXi spans, // spec_pgram
                                  size_t n_groups);




// well_solutions
// Rcpp::NumericVector grf_time_parallel(double radius,
//                                       double storativity,
//                                       double K,
//                                       double thickness,
//                                       const Rcpp::NumericVector& time,
//                                       const Rcpp::NumericVector& flow_rate,
//                                       int flow_time_interval,
//                                       double flow_dimension);
//
// Rcpp::NumericVector hantush_time_parallel(double radius,
//                                           double storativity,
//                                           double transmissivity,
//                                           double leakage,
//                                           const Rcpp::NumericVector &time,
//                                           const Rcpp::NumericVector &flow_rate,
//                                           int flow_time_interval,
//                                           int n_terms);
//
// Rcpp::NumericVector hantush(double radius,
//                             double storativity,
//                             double transmissivity,
//                             double leakage,
//                             const Rcpp::NumericVector &time,
//                             const Rcpp::NumericVector &flow_rate,
//                             const Rcpp::NumericVector &flow_rate_times,
//                             int n_terms);
//
// Rcpp::NumericVector grf(double radius,
//                         double storativity,
//                         double K,
//                         double thickness,
//                         const Rcpp::NumericVector &time,
//                         const Rcpp::NumericVector &flow_rate,
//                         const Rcpp::NumericVector &flow_rate_times,
//                         double flow_dimension);


//==============================================================================
// special_functions.cpp
std::vector<double> bessel_k(std::vector<double> x, int nu);
std::vector<double> bessel_i(std::vector<double> x, int nu);

// Complex bessel functions
Rcpp::ComplexMatrix bessel_i_wrapper(const Rcpp::ComplexVector &x,
                                     double nu,
                                     bool expon_scaled = false,
                                     size_t nseq = 1);
Rcpp::ComplexMatrix bessel_j_wrapper(const Rcpp::ComplexVector &x,
                                     double nu,
                                     bool expon_scaled = false,
                                     size_t nseq = 1);
Rcpp::ComplexMatrix bessel_k_wrapper(const Rcpp::ComplexVector &x,
                                     double nu,
                                     bool expon_scaled = false,
                                     size_t nseq = 1);
Rcpp::ComplexMatrix bessel_y_wrapper(const Rcpp::ComplexVector &x,
                                     double nu,
                                     bool expon_scaled = false,
                                     size_t nseq = 1);

//==============================================================================
// expint.cpp

double expint_E1(double x, int scale);
double expint_E2(double x, int scale);
double gamma_inc(double a, double x);
NumericVector giv(double a, NumericVector x);
NumericVector eiv(NumericVector x, int scale);

NumericVector log_v(NumericVector &x);

//==============================================================================
// well_function_coefficient.cpp
Rcpp::NumericVector well_function_coefficient(const Rcpp::NumericVector flow_rate,
                                              double transmissivity);

double hantush_epsilon(double radius, double leakage);

double grf_coefficient(double flow_rate,
                       double radius,
                       double K,
                       double thickness,
                       double flow_dimension);
Eigen::VectorXd grf_coefficient_eigen(const Eigen::VectorXd &flow_rate,
                                      double radius,
                                      double K,
                                      double thickness,
                                      double flow_dimension);

double grf_u(const double radius,
             const double storativity,
             const double K,
             const double time);

    Rcpp::NumericVector theis_u_time(double radius,
                                     double storativity,
                                     double transmissivity,
                                     const Rcpp::NumericVector &time);

//==============================================================================
// impulse_function.cpp
    std::vector<double> impulse_function(std::vector<double> u,
                                         int flow_time_interval);
