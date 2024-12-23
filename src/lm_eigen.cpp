#include "hydrorecipes.h"


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
Rcpp::List llt_solve_full(Eigen::Map<Eigen::MatrixXd> &X,
                          Eigen::Map<Eigen::MatrixXd> &Y,
                          Rcpp::List subs) {
  const int n(X.rows());
  const int p(X.cols());
  const int n_terms = subs.size();

  const Eigen::LLT<Eigen::MatrixXd> llt(Eigen::MatrixXd(p, p).setZero().selfadjointView<Lower>().
                                          rankUpdate(X.adjoint()));

  const MatrixXd betahat(llt.solve(X.adjoint() * Y));


  const Eigen::MatrixXd fitted(X * betahat);

  Rcpp::List decomposition(n_terms);

  Eigen::VectorXi ind;

  // for (size_t i = 0; i < n_terms; ++i) {
  //   ind = subs[i];
  //   decomposition[i] = X.middleCols(ind[0], ind[1]) * betahat.middleRows(ind[0], ind[1]);
  // }


  const Eigen::MatrixXd resid(Y.array() - fitted.array());

  const unsigned int df(n - p);
  const Eigen::VectorXd s(resid.colwise().norm().array() / std::sqrt(double(df)));
  Eigen::MatrixXd se(s * llt.matrixL().solve(MatrixXd::Identity(p, p)).colwise().norm());

  return Rcpp::List::create(Named("coefficients") = betahat,
                            Named("fitted.values") = fitted,
                            Named("decomposition") = decomposition,
                            Named("residuals") = resid,
                            Named("s") = s,
                            Named("df.residual") = df,
                            Named("rank") = p,
                            Named("Std. Error") = se
  );

}




// [[Rcpp::export]]
Eigen::MatrixXd llt_weighted_solve(Eigen::Map<Eigen::MatrixXd> &X,
                          Eigen::Map<Eigen::MatrixXd> &Y,
                          Eigen::Map<Eigen::VectorXd> &w) {
  const int n(X.rows());
  const int p(X.cols());

  Eigen::MatrixXd out = (X.transpose() * w.asDiagonal() * X).llt().solve(X.transpose() * w.asDiagonal() * Y);

  return(out);
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

x <- matrix(rnorm(10000000), ncol = 50)
y <- matrix(rnorm(2000000), ncol = 10)
w <- rep(1, 200000)
tmp <- hydrorecipes:::llt_weighted_solve(x,y,w)
# yv <- as.numeric(y)

bench::mark(
  tmp <- hydrorecipes:::llt_solve(x, y),
  tmp1 <- hydrorecipes:::llt_fitted(x, y),
  tmpw <- hydrorecipes:::llt_weighted_solve(x, y, w),
  tmp2 <- lm(y~x-1, model = FALSE),
  tmp3 <- lm.fit(x, y),
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
