// [[Rcpp::depends(BH)]]

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppThread)]]

// #define ARMA_DONT_PRINT_ERRORS
// #define ARMA_USE_TBB_ALLOC
//
#define RCPP_ARMADILLO_RETURN_VEC_AS_VECTOR
#define RCPP_ARMADILLO_RETURN_COLVEC_AS_VECTOR
#define RCPP_ARMADILLO_RETURN_ROWVEC_AS_VECTOR

#define EIGEN_FFTW_DEFAULT
#include <unsupported/Eigen/FFT>
#include <unsupported/Eigen/SpecialFunctions>

// #include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/expint.hpp>
#include <boost/math/special_functions/gamma.hpp>
// #include <boost/math/special_functions/factorials.hpp>
// #include <boost/math/special_functions/bessel.hpp>
// #include <boost/math/quadrature/gauss_kronrod.hpp>

#include <Eigen/Eigenvalues>
#include <fftw3.h>
#include <splines2Armadillo.h>
// #include <specialfunctions.h>

#include <RcppEigen.h>
#include <RcppThread.h>


using namespace Rcpp;
using namespace Eigen;
using Eigen::ArrayXd;
using Eigen::ArrayXi;
using Eigen::MatrixXcd;
using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::RowVectorXd;
using Eigen::Vector2d;
using Eigen::Vector3d;
using Eigen::VectorXcd;
using Eigen::VectorXd;
using Eigen::VectorXi;

using Eigen::FFT;


//==============================================================================
// lm_eigen.cpp
Eigen::MatrixXd llt_solve(Eigen::Map<Eigen::MatrixXd> &X,
                          Eigen::Map<Eigen::MatrixXd> &Y);
Eigen::MatrixXd llt_fitted(Eigen::Map<Eigen::MatrixXd> &X,
                           Eigen::Map<Eigen::MatrixXd> &Y);

//==============================================================================
// lags.cpp
int check_lag(int n,
              int lag,
              int n_shift);
int get_length(int n,
               int n_subset);
int get_start(int n_out,
              int lag,
              int n_subset);
int get_end(int n,
            int n_out,
            int lag,
            int n_subset);
NumericVector shift_subset(const NumericVector &x,
                           size_t lag,
                           size_t n_subset,
                           size_t n_shift);
List lag_list(const NumericVector &x,
              const IntegerVector &lags,
              size_t n_subset,
              size_t n_shift);

//==============================================================================
// fft_helpers.cpp
size_t index_from_i_j(size_t i, size_t j, size_t n_col);
size_t index_from_j_i(size_t i, size_t j, size_t n_col);

size_t get_column_number(size_t n);

size_t next_n_eigen(size_t n);
Eigen::VectorXd pad_vector(Eigen::VectorXd x, size_t n_old, size_t n_new);

Eigen::MatrixXd detrend_matrix(const Eigen::MatrixXd &x);
Eigen::MatrixXd demean_matrix(const Eigen::MatrixXd &x);
Eigen::VectorXd detrend_vector(Eigen::VectorXd x);
Eigen::VectorXd demean_vector(Eigen::VectorXd x);
Eigen::MatrixXd detrend_and_demean_matrix(const Eigen::MatrixXd &x,
                                          bool detrend,
                                          bool demean);
Rcpp::List detrend_and_demean_list(Rcpp::List &x,
                                   bool detrend,
                                   bool demean);

Eigen::VectorXd modified_daniell(Eigen::VectorXi spans);
Eigen::MatrixXcd kernel_apply(Eigen::MatrixXcd &x,
                              Eigen::VectorXd &y);
Rcpp::List kernel_apply_list(Rcpp::List x,
                             Eigen::VectorXd &y);
Eigen::VectorXd spec_taper(size_t n_row, double p = 0.1);

Eigen::VectorXi make_groups(size_t n_groups,
                            size_t n);
Eigen::ArrayXd power_spaced(size_t n, double min, double max, double power);
Eigen::VectorXd group_frequency(Eigen::ArrayXd frequencies,
                                size_t n_groups);
Eigen::VectorXd determine_frequency(size_t n);
Eigen::MatrixXcd check_ffts(Eigen::MatrixXcd &x,
                            double cutoff);
Eigen::VectorXi which_indices(const Eigen::VectorXd &x,
                              const Eigen::VectorXd &knots);
Eigen::ArrayXd gamma_inc(Eigen::ArrayXd u, double a);
//==============================================================================
// fft_windows.cpp
Eigen::VectorXd window_hann(size_t n);
Eigen::VectorXcd window_hann_cplx(size_t n);
Eigen::VectorXd window_rectangle(size_t n);
Eigen::VectorXd window_tukey(size_t n, double r);
double window_scale(Eigen::VectorXd window, size_t n_new, size_t n_fft);

//==============================================================================
// fft.cpp
Eigen::MatrixXcd fft_matrix(Eigen::MatrixXd x,
                            size_t n_new);

Eigen::MatrixXcd multiply_ffts(Eigen::MatrixXcd &x);

// convolve
Eigen::VectorXd convolve_vec(Eigen::VectorXd x,
                             Eigen::VectorXd y);
Eigen::VectorXd convolve_overlap_add(Eigen::VectorXd x,
                                     Eigen::VectorXd y);
Eigen::VectorXd convolve_overlap_save(Eigen::VectorXd x,
                                      Eigen::VectorXd y,
                                      int align);
Rcpp::List convolve_overlap_save_list(Eigen::VectorXd x,
                                      Rcpp::List y,
                                      int align);
Eigen::VectorXd convolve_tf(Eigen::VectorXd x,
                            Eigen::VectorXcd y);
Eigen::MatrixXd convolve_matrix(const Eigen::VectorXd &x,
                                const Eigen::MatrixXd &y,
                                bool remove_partial = true,
                                bool reverse = true);
Eigen::VectorXd convolve_filter(const Eigen::VectorXd &x,
                                const Eigen::VectorXd &y,
                                bool remove_partial,
                                bool reverse);
Rcpp::List convolve_list(const Eigen::VectorXd &x,
                         const List y,
                         const bool remove_partial,
                         const bool reverse);
std::list<Eigen::VectorXd> convolve_list2(const Eigen::VectorXd &x,
                                          const std::list<Eigen::VectorXd> y,
                                          const bool remove_partial,
                                          const bool reverse);
// Spectrum
Eigen::MatrixXcd spec_welch(Eigen::MatrixXd &x,
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

Eigen::MatrixXcd spec_pgram(Eigen::MatrixXd &x,
                            const Eigen::VectorXi &spans,
                            bool detrend,
                            bool demean,
                            double taper = 0.1,
                            bool pad_fft = true);

// Solve

Eigen::MatrixXcd solve_cplx_parallel(const Eigen::MatrixXcd &x);
Eigen::MatrixXcd solve_cplx_irr(Eigen::MatrixXcd &x,
                                size_t n_groups);

Eigen::MatrixXcd transfer_pgram_smooth(Eigen::MatrixXd &x,
                                       const Eigen::VectorXi &spans,
                                       bool detrend,
                                       bool demean,
                                       double taper,
                                       size_t n_groups);
Eigen::MatrixXcd transfer_pgram(Eigen::MatrixXd &x,
                                const Eigen::VectorXi &spans,
                                bool detrend,
                                bool demean,
                                double taper,
                                size_t n_groups);

Eigen::MatrixXcd transfer_welch(Eigen::MatrixXd &x,
                                size_t length_subset,
                                double overlap,
                                Eigen::VectorXd window = Eigen::VectorXd::Zero(0));

// Processing
Eigen::MatrixXd ordinary_coherence_phase(const Eigen::MatrixXcd &x);
Eigen::MatrixXd frequency_to_time_domain(Eigen::MatrixXcd &pgram,
                                         size_t n_groups);
// Eigen::ArrayXd frf_to_brf(const Eigen::VectorXcd& x,
//                           std::complex<double> dc1,
//                           std::complex<double> dc2);
// Eigen::MatrixXcd interpolate_tf(Eigen::MatrixXcd& x,
//                                  const Eigen::ArrayXd& frequency_irregular,
//                                  const Eigen::ArrayXd& frequency_regular,
//                                  Eigen::VectorXd& knots);
// Eigen::MatrixXd predict_pgram_frf(Eigen::MatrixXd& x,
//                                   Eigen::VectorXi span);
Eigen::MatrixXd predict_pgram_frf(Eigen::MatrixXd &x,
                                  Eigen::MatrixXd &x_out,
                                  Eigen::VectorXi spans, // spec_pgram
                                  size_t n_groups);

//==============================================================================
// b_spline_arma.cpp
Rcpp::List b_spline_list(const arma::vec &x,
                         const unsigned int df,
                         const unsigned int degree,
                         const arma::vec &internal_knots,
                         const arma::vec &boundary_knots,
                         const bool complete_basis,
                         const bool periodic,
                         const unsigned int derivs,
                         const bool integral);
Rcpp::List n_spline_list(const arma::vec &x,
                         const unsigned int df,
                         const unsigned int degree,
                         const arma::vec &internal_knots,
                         const arma::vec &boundary_knots,
                         const bool complete_basis,
                         const bool periodic,
                         const unsigned int derivs,
                         const bool integral);
Rcpp::List b_spline_list2(const arma::vec &x,
                          const unsigned int df,
                          const unsigned int degree,
                          const arma::vec &internal_knots,
                          const arma::vec &boundary_knots,
                          const bool complete_basis,
                          const bool periodic,
                          const unsigned int derivs,
                          const bool integral);
std::list<Eigen::VectorXd> b_spline_list3(const arma::vec &x,
                                          const unsigned int df,
                                          const unsigned int degree,
                                          const arma::vec &internal_knots,
                                          const arma::vec &boundary_knots,
                                          const bool complete_basis,
                                          const bool periodic,
                                          const unsigned int derivs,
                                          const bool integral);
arma::vec log_lags_arma(arma::uword n, arma::uword lag_max);

//==============================================================================
// distributed_lag.cpp
Eigen::MatrixXd distributed_lag_thread(const Eigen::VectorXd &x,
                                       const Eigen::MatrixXd &bl,
                                       int n_thread);
Rcpp::List distributed_lag_eigen(Eigen::Map<Eigen::VectorXd> x,
                                 Eigen::Map<Eigen::MatrixXd> bl);
List distributed_lag_thread_eigen(Eigen::Map<Eigen::VectorXd> x,
                                  Eigen::Map<Eigen::MatrixXd> bl,
                                  int lag_max,
                                  int n_subset,
                                  int n_shift,
                                  int n_thread);
Rcpp::List distributed_lag_list(Eigen::Map<Eigen::VectorXd> x,
                                arma::uword n_lag,
                                arma::uword lag_max,
                                const unsigned int df,
                                const unsigned int degree,
                                const arma::vec &internal_knots,
                                const arma::vec &boundary_knots,
                                const bool complete_basis = true,
                                const bool periodic = false,
                                const unsigned int derivs = 0,
                                const bool integral = false);
Rcpp::List distributed_lag_list3(Eigen::VectorXd x,
                                 arma::uword n_lag,
                                 arma::uword lag_max,
                                 const unsigned int df,
                                 const unsigned int degree,
                                 const arma::vec &internal_knots,
                                 const arma::vec &boundary_knots,
                                 const bool complete_basis = true,
                                 const bool periodic = false,
                                 const unsigned int derivs = 0,
                                 const bool integral = false);

//==============================================================================
// helpers_aquifer.cpp
double dimensionless_well_bore_storage(double rc,
                                       double rw,
                                       double Ss,
                                       double l,
                                       double d);


//==============================================================================
// laplace_classes.cpp
Rcpp::ComplexMatrix bessel_k_cplx(const Rcpp::ComplexMatrix &x,
                                  double nu,
                                  bool expon_scaled,
                                  size_t n_seq);
