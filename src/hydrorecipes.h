// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppThread)]]

#include <RcppEigen.h>
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

using namespace Rcpp;


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
Eigen::MatrixXd distributed_lag_parallel(const Eigen::VectorXd& x,
                             const Eigen::MatrixXd& bl,
                             int lag_max,
                             int n_subset = 1,
                             int n_shift = 0);
