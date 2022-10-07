#define ARMA_DONT_PRINT_ERRORS
// #define ARMA_USE_TBB_ALLOC

#include <RcppArmadillo.h>
#include <RcppParallel.h>

using namespace RcppParallel;
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// to_dummy.cpp
IntegerMatrix to_dummy(const IntegerVector& x,
                       size_t n_fact);


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
mat distributed_lag_parallel(const vec& x,
                             const mat& bl,
                             int lag_max,
                             int n_subset = 1,
                             int n_shift = 0);
mat arma_shift(mat x, int n);
