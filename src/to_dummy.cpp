#include "frecipes.h"


// [[Rcpp::export]]
IntegerVector fi(const NumericVector& x,
                 const NumericVector& vec,
                 const bool rightmost_closed,
                 const bool all_inside,
                 const bool left_open) {

  Rcpp::Function f("findInterval");

  const IntegerVector out = f(x, vec, rightmost_closed, all_inside, left_open);

  return(out);

}


// [[Rcpp::export]]
Rcpp::IntegerVector to_dummy_list_base(const Rcpp::IntegerVector& x,
                                       const int n_fact) {

  Rcpp::IntegerVector z = Rcpp::clone(x);
  z.attr("levels") = R_NilValue;
  z.attr("class") = R_NilValue;
  for (auto& elem: z) {
    if (elem == n_fact) {
      elem = 1;
    } else {
      elem = 0;
    }
  }

  return z;
}


//==============================================================================
//' @title
//' to_dummy
//'
//' @description
//' Create binary terms based on a factor column.
//'
//' @params ind integer vector of values to dummy encode
//'
//' @return List of dummy encoded terms
//'
//' @export
//'
//'
// [[Rcpp::export]]
List to_dummy(const IntegerVector& ind, const bool one_hot) {

  List out;
  const IntegerVector fact = sort_unique(ind);
  const size_t n = fact.size();

  size_t start = 0;
  if (!one_hot & (n > 1)) {
    start = 1;
  }

  for (size_t i = start; i < n; ++i) {
    out.push_back(to_dummy_list_base(ind, fact[i]), std::to_string(fact[i]));
  }

  return(out);

}

//==============================================================================
//' @title
//' to_dummy_list
//'
//' @description
//' Create binary terms based on intervals. This function uses `findInterval`,
//' followed by a conversion to dummy encoding.
//'
//' @inheritParams findInterval
//'
//' @return List of dummy encoded terms
//'
//' @export
//'
//'
// [[Rcpp::export]]
List to_dummy_list(const NumericVector& x,
                   const NumericVector& vec,
                   const bool one_hot = false,
                   const bool rightmost_closed = false,
                   const bool all_inside = false,
                   const bool left_open = false
                   ) {

  const IntegerVector ind = fi(x,
                               vec,
                               rightmost_closed,
                               all_inside,
                               left_open);


  return(to_dummy(ind, one_hot));

}




// // [[Rcpp::export]]
// IntegerMatrix to_dummy7(const IntegerVector& x,
//                         size_t n_fact,
//                         bool intercept = true) {
//
//   IntegerMatrix y(x.size(), n_fact);
//
//   if (!intercept) {
//     n_fact -= 1;
//   }
//
//   for (size_t i = 0; i < n_fact; ++i) {
//     y.column(i) = x == i;
//   }
//
//   return y;
// }
//
// // [[Rcpp::export]]
// IntegerMatrix to_dummy2(const IntegerVector& x,
//                         int n_fact,
//                         bool intercept = true) {
//
//   IntegerMatrix y(x.size(), n_fact);
//
//   if (!intercept) {
//     n_fact -= 1;
//   }
//
//   for (int i = 0; i < x.size(); ++i) {
//     y(i, x(i)) = 1;
//   }
//
//   return y;
// }
//
// // [[Rcpp::export]]
// IntegerVector to_dummy0(const IntegerVector& x,
//                         size_t n_fact,
//                         bool intercept = true) {
//
//   int n = x.size();
//   // defaults to zero
//   IntegerMatrix y(x.size(), n_fact);
//   IntegerVector z = x * n + Rcpp::Range(0, n-1);
//   // one traverse of the vector
//   for (auto& elem: z) {
//
//     y[elem] = 1;
//
//   }
//
//   return x;
// }
//
// // [[Rcpp::export]]
// IntegerMatrix to_dummy(const IntegerVector& x,
//                        size_t n_fact) {
//
//   // defaults to zero
//   IntegerMatrix y(x.size(), n_fact);
//
//   auto i = 0;
//
//   // one traverse of the vector
//   for (auto& elem: x) {
//
//     y(i, elem - 1) = 1;
//     i += 1;
//
//   }
//
//   return y;
// }
//
// // [[Rcpp::export]]
// arma::sp_mat to_dummy13(arma::urowvec& x,
//                         int n_fact) {
//   int n = x.n_elem;
//   arma::umat inds(2, n);
//
//   inds.row(1) = x;
//   inds.row(0) = arma::regspace<arma::urowvec>(0, n - 1);
//
//   arma::vec z(n);
//   z.ones();
//
//   arma::sp_mat y(inds, z);
//
//   return(y);
//
// }
//
// // [[Rcpp::export]]
// arma::imat to_dummy3(arma::ivec& x,
//                      int n_fact) {
//
//   int n = x.n_elem;
//   arma::imat y(n, n_fact);
//
//   for (int i = 0; i < n_fact; ++i) {
//     y.elem(arma::find(x == i)).ones();
//   }
//
//   return y;
// }
//
// // [[Rcpp::export]]
// imat to_dummy4(const ivec& x,
//                int n_fact,
//                bool intercept = true) {
//
//   int n = x.n_elem;
//   imat y(n, n_fact);
//
//   if (!intercept) {
//     n_fact -= 1;
//   }
//
//   for (int i = 0; i < n; ++i) {
//     y(i, x(i)) = 1;
//   }
//
//   return y;
// }
//
// // [[Rcpp::export]]
// imat to_dummy10(const ivec& x,
//                 int n_fact,
//                 bool intercept = true) {
//
//   int n = x.n_elem;
//   imat y(n_fact,n);
//
//   if (!intercept) {
//     n_fact -= 1;
//   }
//
//   for (int i = 0; i < n; ++i) {
//     y(x(i), i) = 1;
//   }
//
//   return y;
// }
//
// // [[Rcpp::export]]
// imat to_dummy5(const arma::ivec& x,
//                size_t n_fact,
//                bool intercept = true) {
//
//   size_t n = x.n_elem;
//   arma::imat y( n_fact, n);
//
//   if (!intercept) {
//     n_fact -= 1;
//   }
//
//   for (size_t i = 0; i < n; ++i) {
//     y(x(i), i) = 1;
//   }
//
//   return y.t();
// }
//
//
// // // [[Rcpp::export]]
// // imat to_dummy6(const arma::ivec& x,
// //                size_t n_fact,
// //                bool intercept = true) {
// //
// //   size_t n = x.n_elem;
// //   arma::imat y( n, n_fact);
// //   y.fill(0);
// //   arma::ivec z(n);
// //   if (!intercept) {
// //     n_fact -= 1;
// //   }
// //
// //   for (uword i = 0; i < n_fact; ++i) {
// //     y(find(x==i)).ones();
// //   }
// //
// //   return y;
// // }
//
// // [[Rcpp::export]]
// Eigen::MatrixXi to_dummy8(Eigen::VectorXi x,
//                           size_t n_fact,
//                           bool intercept = true) {
//
//   int n = x.size();
//   Eigen::MatrixXi y(n, n_fact);
//
//   if (!intercept) {
//     n_fact -= 1;
//   }
//
//   for (size_t i = 0; i < n; ++i) {
//     y(i, x(i)) = 1;
//   }
//
//   return y;
// }
//
// // [[Rcpp::export]]
// Eigen::MatrixXi to_dummy9(const Eigen::VectorXi& x,
//                           size_t n_fact,
//                           bool intercept = true) {
//
//   int n = x.size();
//   Eigen::MatrixXi y(n_fact, n);
//
//   if (!intercept) {
//     n_fact -= 1;
//   }
//
//   for (size_t i = 0; i < n; ++i) {
//     y(x(i), i) = 1;
//   }
//
//   return y;
// }
//
// // [[Rcpp::export]]
// Eigen::MatrixXi to_dummy11(Eigen::VectorXi x,
//                            size_t n_fact,
//                            bool intercept = true) {
//
//   int n = x.size();
//   int ind;
//
//   Eigen::MatrixXi y(n, n_fact);
//
//   if (!intercept) {
//     n_fact -= 1;
//   }
//   x = n*x;
//   for (int i = 0; i < n; ++i) {
//     ind = i + x(i);
//     y(ind) = 1;
//   }
//
//   return y;
// }
//
// // [[Rcpp::export]]
// List to_dummy_lst(const IntegerVector& x,
//                   size_t n_fact,
//                   CharacterVector nms,
//                   bool intercept = true) {
//
//   List out_lst;
//   IntegerVector y(x.length());
//   std::string nm;
//
//   if (!intercept) {
//     n_fact -= 1;
//   }
//
//   for (size_t i = 0; i < n_fact; i++) {
//     y = (x == i);
//     nm = nms(i);
//     out_lst.push_back(clone(y), nm);
//   }
//
//   return out_lst;
// }
//
//
//
// // [[Rcpp::export]]
// List to_dummy_lst2(const IntegerMatrix& x,
//                    size_t n_fact,
//                    CharacterVector nms,
//                    bool intercept = true) {
//
//   IntegerVector y(x.nrow());
//
//   std::string nm;
//
//   if (!intercept) {
//     n_fact -= 1;
//   }
//
//   List out_lst;
//   for (size_t j = 0; j < x.ncol(); j++) {
//     for (size_t i = 0; i < n_fact; i++) {
//       y = (x.column(j) == i);
//       nm = nms(i);
//       out_lst.push_back(clone(y), nm);
//     }
//   }
//
//   return out_lst;
// }
//
// // [[Rcpp::export]]
// List to_dummy_lst3(arma::irowvec& x,
//                    const size_t n_fact) {
//
//   size_t n = x.size();
//   arma::irowvec y(n);
//
//   x -= 1;
//   List out_lst;
//   for (size_t i = 0; i < n_fact; i++) {
//     y.zeros();
//     y.elem(arma::find(x == i)).ones();
//
//     out_lst.push_back(y);
//   }
//
//   return out_lst;
// }
//
// // [[Rcpp::export]]
// List to_dummy_lst4(arma::irowvec& x,
//                    const size_t n_fact) {
//
//   size_t n = x.size();
//
//   // List out_lst(n_fact);
//   arma::irowvec y(n);
//   std::list<arma::irowvec> test;
//
//   for (size_t i = 0; i < n_fact; i++) {
//     y.zeros();
//     y.elem(arma::find(x == i)).ones();
//     test.push_back(y);
//   }
//   Rcpp::List result = Rcpp::List::create(test);
//   return result;
// }
//
//
// // [[Rcpp::export]]
// List to_dummy_lst5(const Rcpp::IntegerVector& x,
//                          const int n_fact) {
//
//   List out;
//
//   for (int i = 1; i <= n_fact; ++i) {
//     out.push_back(to_dummy_lst_base(x, i));
//   }
//
//   return out;
// }
//
// // [[Rcpp::export]]
// List to_dummy_df(const IntegerVector& x,
//                  size_t n_fact,
//                  bool intercept = true) {
//
//   // IntegerVector y(x.length());
//
//   if (!intercept) {
//     n_fact -= 1;
//   }
//   List out_lst(n_fact);
//
//   for (size_t i = 0; i < n_fact; i++) {
//     // y = (x == i);
//     out_lst(i) = (x == i);
//
//   }
//
//   return (out_lst);
// }





/*** R
set.seed(123)
x <- sort(rnorm(2e7))
vec = -7:7

bench::mark(
  frecipes:::to_dummy_list(x, vec),
  iterations = 10
)



*/
