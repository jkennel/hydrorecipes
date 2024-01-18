#include "frecipes.h"



// [[Rcpp::export]]
Eigen::MatrixXd llt_solve(Eigen::Map<Eigen::MatrixXd> &X,
                          Eigen::Map<Eigen::MatrixXd> &Y) {
  const int n(X.rows());
  const int p(X.cols());

  const Eigen::LLT<Eigen::MatrixXd> llt(Eigen::MatrixXd(p, p).setZero().selfadjointView<Lower>().
                                          rankUpdate(X.adjoint()));
  const MatrixXd betahat(llt.solve(X.adjoint() * Y));
  // const Eigen::MatrixXd fitted(X * betahat);
  // const Eigen::MatrixXd resid(Y.array() - fitted.array());
  //
  // const unsigned int df(n - p);
  // const Eigen::VectorXd s(resid.colwise().norm().array() / std::sqrt(double(df)));
  // Eigen::MatrixXd se(s * llt.matrixL().solve(MatrixXd::Identity(p, p)).colwise().norm());
  // return Rcpp::List::create(Named("coefficients") = betahat,
  //                           Named("fitted.values") = fitted,
  //                           Named("residuals") = resid,
  //                           Named("s") = s,
  //                           Named("df.residual") = df,
  //                           Named("rank") = p,
  //                           Named("Std. Error") = se
  // );

  return(betahat);
}

// [[Rcpp::export]]
Eigen::MatrixXd llt_fitted(Eigen::Map<Eigen::MatrixXd> &X,
                          Eigen::Map<Eigen::MatrixXd> &Y) {
  const int n(X.rows());
  const int p(X.cols());

  const Eigen::LLT<Eigen::MatrixXd> llt(Eigen::MatrixXd(p, p).setZero().selfadjointView<Lower>().
                                          rankUpdate(X.adjoint()));
  const MatrixXd betahat(llt.solve(X.adjoint() * Y));
  const Eigen::MatrixXd fitted(X * betahat);
  // const Eigen::MatrixXd resid(Y.array() - fitted.array());
  //
  // const unsigned int df(n - p);
  // const Eigen::VectorXd s(resid.colwise().norm().array() / std::sqrt(double(df)));
  // Eigen::MatrixXd se(s * llt.matrixL().solve(MatrixXd::Identity(p, p)).colwise().norm());
  // return Rcpp::List::create(Named("coefficients") = betahat,
  //                           Named("fitted.values") = fitted,
  //                           Named("residuals") = resid,
  //                           Named("s") = s,
  //                           Named("df.residual") = df,
  //                           Named("rank") = p,
  //                           Named("Std. Error") = se
  // );

  return(fitted);
}

/*** R

x <- matrix(rnorm(100000000), ncol = 50)
y <- matrix(rnorm(20000000), ncol = 10)
# yv <- as.numeric(y)

bench::mark(
  tmp <- frecipes:::llt_solve(x,y),
  tmp2 <- lm(y~x-1, model = FALSE),
  a <- x[,2:3] %*% tmp[2:3,],
  check = FALSE
)

s <- c(1, 10, 20, 30, 40)
e <- c(9, 19, 29, 39, 49)

bench::mark(
  for ( i in seq_along(s)) {
    a <- x[,s[i]:e[i]] %*% tmp[s[i]:e[i],]
  }
)

*/
