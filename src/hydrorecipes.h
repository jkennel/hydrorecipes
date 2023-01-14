// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppThread)]]

#define EIGEN_FFTW_DEFAULT

#include <RcppEigen.h>
#include <unsupported/Eigen/FFT>
#include <unsupported/Eigen/Splines>
#include <fftw3.h>
#include <Eigen/StdVector>

#include <RcppThread.h>

using namespace Eigen;
using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::ArrayXd;
using Eigen::ArrayXi;
using Eigen::MatrixXcd;
using Eigen::VectorXcd;
using Eigen::Vector3d;
using Eigen::Vector2d;
using Eigen::RowVectorXd;

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

Eigen::VectorXd modified_daniell(Eigen::VectorXi spans, size_t n);
Eigen::MatrixXcd kernel_apply(Eigen::MatrixXcd& x,
                              Eigen::VectorXd& y);
Eigen::VectorXd spec_taper(size_t n_row, double p = 0.1);



Eigen::VectorXi make_groups(size_t n_row,
                            size_t power,
                            size_t n_groups,
                            size_t min_aggregate);
Eigen::ArrayXd power_spaced(size_t n, double min, double max, double power);
Eigen::VectorXd group_frequency(Eigen::ArrayXd frequencies,
                                size_t power,
                                size_t n_groups,
                                size_t min_aggregate);
Eigen::VectorXd determine_frequency(size_t n);
Eigen::MatrixXcd check_ffts(Eigen::MatrixXcd& x,
                            double cutoff);
Eigen::VectorXi which_indices(const Eigen::VectorXd& x,
                              const Eigen::VectorXd& knots);



//==============================================================================
// b_spline.cpp
Eigen::MatrixXd b_spline(const Eigen::ArrayXd& x,
                         const Eigen::ArrayXd& knots,
                         size_t degree);


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
Eigen::VectorXd convolve_tf(Eigen::VectorXd x, Eigen::VectorXcd y);
Eigen::MatrixXd convolve_matrix(const Eigen::VectorXd& x,
                                const Eigen::MatrixXd& y,
                                bool remove_partial = true,
                                bool reverse = true);



// Spectrum
Eigen::MatrixXcd spec_welch(const Eigen::MatrixXd& x,
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
                                size_t power,
                                size_t n_groups,
                                size_t min_aggregate);

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
Eigen::MatrixXd frequency_to_time_domain(Eigen::MatrixXcd& x);
// Eigen::ArrayXd frf_to_brf(const Eigen::VectorXcd& x,
//                           std::complex<double> dc1,
//                           std::complex<double> dc2);
Eigen::MatrixXcd interpolate_frf(const Eigen::ArrayXd& x,
                                 const Eigen::MatrixXcd& y,
                                 Eigen::VectorXd& knots,
                                 size_t degree,
                                 const Eigen::ArrayXd& x_interp);
// Eigen::MatrixXd predict_pgram_frf(Eigen::MatrixXd& x,
//                                   Eigen::VectorXi span);
// Eigen::MatrixXd predict_pgram_smooth_frf(Eigen::MatrixXd& x,
//                                          Eigen::VectorXi span,
//                                          double power,
//                                          size_t n_groups);
