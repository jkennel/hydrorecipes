#include "frecipes.h"


// [[Rcpp::export]]
int check_lag(int n,
              int lag,
              int n_shift) {

  if (lag - n_shift > n) {
    throw std::range_error(
        "lag + n_shift cannot be greater than series length");
  }
  if (n_shift == 0) {
    return(lag);
  } else {
    return(lag - n_shift);
  }

}

// [[Rcpp::export]]
int get_length(int n,
               int n_subset) {

  int n_out;

  // Get length
  if (n_subset == 1){
    n_out = n;
  } else {
    if ((n % n_subset) != 0) {
      n_out = (n / n_subset) + 1;
    } else {
      n_out = (n / n_subset);
    }
  }
  return(n_out);
}


// [[Rcpp::export]]
int get_start(int n_out,
              int lag,
              int n_subset) {

  int start;

  // bounds
  if ((lag % n_subset) != 0) {
    if (lag > 0) {
      start = (lag / n_subset) + 1;
    } else {
      start = 0;
    }
  } else {
    if (lag > 0) {
      start = (lag / n_subset);
    } else {
      start = 0;
    }
  }

  return(start);
}

// [[Rcpp::export]]
int get_end(int n,
            int n_out,
            int lag,
            int n_subset) {

  int end;
  // bounds
  if ((lag % n_subset) != 0) {
    if (lag > 0) {
      end = n_out;
    } else {
      if ((n % n_subset) != 0) {
        end = n_out - (-lag / n_subset + 1);
      } else {
        end = n_out - (-lag / n_subset);
      }
    }
  } else {
    if (lag > 0) {
      end = n_out;
    } else {
      end = n_out - (-lag / n_subset);
    }
  }
  return(end);
}

//' @title
//' shift_subset
//'
//' @description
//' lag data and subset the results
//'
//' @inheritParams step_lead_lag
//' @param x to lag (numeric vector)
//' @param lag amount to lag or lead if negative (integer)
//'
//' @return vector with lagged values
//'
//' @noRd
//'
// [[Rcpp::export]]
Rcpp::NumericVector shift_subset(const Rcpp::NumericVector& x,
                                 size_t lag,
                                 size_t n_subset,
                                 size_t n_shift) {

  if (n_shift >= n_subset) {
    throw std::range_error("shift_subset: n_shift must be less than n_subset");
  }

  int n = x.size();
  int n_out;
  int start, end;
  int wh;


  lag   = check_lag(n, lag, n_shift);
  n_out = get_length(n, n_subset);

  start = get_start(n_out, lag, n_subset);
  end   = get_end(n, n_out, lag, n_subset);

  Rcpp::NumericVector out(n_out, NA_REAL);

  if (start >= end) {
    throw std::range_error("shift_subset: the number of lags, n_subset or n_shift is too large");
  }

  for (int i = start; i < end; ++i) {
    wh = (i * n_subset) - lag;

    out[i] = x[wh];

  }

  return(out);
}


//==============================================================================
//' @title
//' lag_list
//'
//' @description
//' Create lagged terms
//'
//' @param x numeric vector - variable to lag
//' @param lags integer vector - amount to lag
//' @param n_subset take every n_subset rows
//' @param n_shift shift values from starting on first row.  Should be less than
//'  n_subset
//'
//' @return List of lagged terms
//'
//' @export
//'
//' @noRd
//'
// [[Rcpp::export]]
Rcpp::List lag_list(const Rcpp::NumericVector& x,
              const Rcpp::IntegerVector& lags,
              size_t n_subset,
              size_t n_shift
) {

  size_t n = x.size();
  size_t n_row;

  if (n_subset == 1){
    n_row = (n - n_shift);
  } else {
    n_row = ((n - n_shift - 1) / n_subset) + 1;
  }
  size_t n_col = lags.size();

  Rcpp::List out(n_col);

  for (size_t i = 0; i < n_col; ++i) {
    out[i] = shift_subset(x, lags(i), n_subset, n_shift);
  }


  return(out);
}


// //' @title
//  //' lag_matrix
//  //'
//  //' @description
//  //' lag data and subset the results
//  //'
//  //' @inheritParams step_lead_lag
//  //' @param x to lag (numeric vector)
//  //' @param lags lead or lag values (numeric vector)
//  //' @param var_name name for the generated matrix columns (character)
//  //'
//  //' @return matrix with lagged values
//  //'
//  //' @noRd
//  // [[Rcpp::export]]
//  NumericMatrix lag_matrix(const NumericMatrix& x,
//                           const IntegerVector& lags,
//                           CharacterVector suffix,
//                           std::string prefix,
//                           int n_subset,
//                           int n_shift
//  ) {
//
//    int n = x.nrow();
//    int n_row;
//
//    if (n_subset == 1){
//      n_row = (n - n_shift);
//    } else {
//      n_row = ((n - n_shift - 1) / n_subset) + 1;
//    }
//    int n_col = lags.size();
//    int n_var = x.ncol();
//
//    CharacterVector nm(n_col * n_var);
//    NumericMatrix out = NumericMatrix(n_row, n_col * n_var);
//
//    for (int j = 0; j < n_var; ++j) {
//      for (int i = 0; i < n_col; ++i) {
//        out(_, i + j * n_col) = shift_subset(x(_, j), lags[i], n_subset, n_shift);
//
//        // // Column names
//        // if (lags[i] < 0) {
//        //   nm[i + j * n_col] = prefix + 'n' + std::to_string(abs(lags[i])) + '_' + suffix[j];
//        // } else {
//        //   nm[i + j * n_col] = prefix + std::to_string(lags[i]) + '_' + suffix[j];
//        // }
//
//      }
//    }
//
//    // colnames(out) = nm;
//
//    return(out);
//  }
// // [[Rcpp::export]]
// colvec arma_shift(const colvec& x, int n) {
//
//   if (n == 0) {
//     return(x);
//   }
//
//   colvec out = shift(x, n);
//   colvec fill_vec(abs(n), fill::value(NA_REAL));
//
//   if (n > 0) {
//     out.head(n) = fill_vec;
//   } else {
//     out.tail(abs(n)) = fill_vec;
//   }
//
//   return out;
// }
//
// // [[Rcpp::export]]
// colvec arma_shift_subset(const colvec& x,
//                          size_t lag,
//                          size_t n_subset,
//                          size_t n_shift) {
//
//   if (n_shift >= n_subset) {
//     throw std::range_error("shift_subset: n_shift must be less than n_subset");
//   }
//
//   std::size_t n = x.n_elem;
//   std::size_t n_out;
//   std::size_t start, end;
//   // std::size_t wh;
//
//
//   lag   = check_lag(n, lag, n_shift);
//   n_out = get_length(n, n_subset);
//
//   start = get_start(n_out, lag, n_subset);
//   end   = get_end(n, n_out, lag, n_subset);
//
//   if (start >= end) {
//     throw std::range_error("shift_subset: the number of lags, n_subset or n_shift is too large");
//   }
//
//   colvec out(n_out, fill::value(NA_REAL));
//   uvec inds = regspace<uvec>((start * n_subset) - lag,
//                              n_subset,
//                              (end * n_subset) - lag - 1);
//
//   out.tail(inds.n_elem) = x.elem(inds);
//
//   return(out);
// }
//
// // [[Rcpp::export]]
// field<colvec> lag_lead_matrix(const colvec& x,
//                               ivec n,
//                               size_t n_subset,
//                               size_t n_shift) {
//
//   field<colvec> out(n.n_elem);
//
//   for (size_t i=0; i < n.n_elem; ++i) {
//
//     if (n_subset == 1 & n_shift == 0) {
//       out[i] = arma_shift(x, n[i]);
//     } else {
//       out[i] = arma_shift_subset(x, n[i], n_subset, n_shift);
//     }
//
//   }
//
//   return(out);
// }
//
// // [[Rcpp::export]]
// mat f_u(mat x){
//   return flipud(x);
// }
//
// // [[Rcpp::export]]
// mat f_l(mat x){
//   return fliplr(x);
// }
/*** R

n <- 2000000L
m <- 1:n
nn <- 5
bench::mark(
  m1 <- lag_list(m, -1:nn, n_subset = 1, n_shift = 0),
  check = FALSE,
  iterations = 3
)


a <- list()
for ( i in 1:10) {
  a[[i]] <- data.frame(x = rnorm(1e6))
}
z <- list()
bench::mark(
  q <- lapply(a, FUN=function(x) x*2),
  {for (i in 1:length(a)) { z[[i]] = 2*a[[i]]}},
  check = FALSE
)

microbenchmark::microbenchmark(
  q <- lapply(a, FUN=function(x) x*2),
  {for (i in 1:length(a)) { z[[i]] = 2*a[[i]]}}
)
*/
