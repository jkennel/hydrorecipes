// // [[Rcpp::depends(RcppArmadillo)]]
//
// #include <RcppArmadillo.h>
//
// // [[Rcpp::export]]
// Rcpp::List drop_rows(const Rcpp::List& x,
//                      const Rcpp::IntegerVector& vec) {
//
//   size_t n = x.length();
//   size_t n_sub = vec.size();
//   Rcpp::List out;
//   Rcpp::NumericVector l_vec(n);
//
//   for (size_t i = 0; i < n; ++i) {
//     l_vec = x[i];
//     out.push_back(l_vec[vec]);
//   }
//
//   return(out);
// }
//
// // [[Rcpp::export]]
// arma::vec drop_rows_base(const arma::vec& x,
//                     const arma::uvec& vec) {
//
//   size_t n = x.size();
//   arma::vec out = x.elem(vec);
//
//   return(out);
// }
//
// // [[Rcpp::export]]
// Rcpp::List drop_rows_list(const Rcpp::List& x,
//                          const arma::uvec& vec) {
//
//   size_t n = x.size();
//   Rcpp::List out;
//
//   for (size_t i = 0; i < n; ++i) {
//     // l_vec = x[i];
//     out.push_back(drop_rows_base(x[i], vec));
//   }
//
//   return(out);
// }
//
// /*** R
//
// library(data.table)
// n <- 1e7
// df <- data.frame(age=sample(1:65,1e7,replace=TRUE),x=rnorm(1e7),y=rpois(1e7,25))
// dt <- as.data.table(df)
// vec <- df[,1]
// subsetter <- function(x, y) {
//   lapply(unclass(x), '[', -y)
// }
// to_keep <- seq(1, 1e7, 100)
// bench::mark(
//   subsetter(df, to_keep),
//   drop_rows(unclass(df), to_keep),
//   drop_rows_base(vec, to_keep),
//   vec[to_keep],
//   dt[to_keep],
//   df[to_keep,],
//   qDT(df)[to_keep],
//   check = FALSE
// )
//
// bench::mark(
//   # drop_rows(a$a, to_keep),
//   drop_rows_list(unlist(tmp$result, recursive = FALSE), to_keep),
// )
// */
