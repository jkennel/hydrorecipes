#include "hydrorecipes.h"

// [[Rcpp::export]]
IntegerMatrix to_dummy(IntegerVector x,
                       size_t n_fact) {

  // factor to c++ indexing
  x = x - 1;

  // defaults to zero
  IntegerMatrix y(x.size(), n_fact);

  // one traverse of the vector
  for (auto& elem: x) {

    auto i = &elem - &x[0];
    y(i, elem) = 1;

  }


  return y;
}


/*** R
*/
