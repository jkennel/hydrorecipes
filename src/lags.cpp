// inspired by http://davegiles.blogspot.com/2017/01/explaining-almon-distributed-lag-model.html
// and the dlnm package

#define ARMA_DONT_PRINT_ERRORS
// #define ARMA_USE_TBB_ALLOC

#include <RcppArmadillo.h>
#include <RcppParallel.h>

using namespace RcppParallel;
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
int check_lag(int n, int lag, int n_shift) {

  if(lag - n_shift > n) {
    throw std::range_error(
        "lag + n_shift cannot be greater than series length");
  }
  if(n_shift == 0) {
    return(lag);
  } else {
    return(lag - n_shift);
  }

}

// [[Rcpp::export]]
int get_length(int n, int n_subset) {

  int n_out;

  // Get length
  if(n_subset == 1){
    n_out = n;
  } else {
    if((n % n_subset) != 0) {
      n_out = (n / n_subset) + 1;
    } else {
      n_out = (n / n_subset);
    }
  }
  return(n_out);
}


// [[Rcpp::export]]
int get_start(int n_out, int lag, int n_subset) {
  int start;

  // bounds
  if((lag % n_subset) != 0) {
    if(lag > 0) {
      start = (lag / n_subset) + 1;
    } else {
      start = 0;
    }
  } else {
    if(lag > 0) {
      start = (lag / n_subset);
    } else {
      start = 0;
    }
  }

  return(start);
}

// [[Rcpp::export]]
int get_end(int n, int n_out, int lag, int n_subset) {

  int end;
  // bounds
  if((lag % n_subset) != 0) {
    if(lag > 0) {
      end = n_out;
    } else {
      if((n % n_subset) != 0) {
        end = n_out - (-lag / n_subset + 1);
      } else {
        end = n_out - (-lag / n_subset);
      }
    }
  } else {
    if(lag > 0) {
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
Rcpp::NumericVector shift_subset(const Rcpp::NumericVector x,
                                 int lag = 0,
                                 int n_subset = 1,
                                 int n_shift = 0) {

  int n = x.size();
  int n_out;
  int start, end;
  int wh;


  lag   = check_lag(n, lag, n_shift);
  n_out = get_length(n, n_subset);

  Rcpp::NumericVector out(n_out, NA_REAL);

  start = get_start(n_out, lag, n_subset);
  end   = get_end(n, n_out, lag, n_subset);

  for (int i = start; i < end; i++) {
    wh = (i * n_subset) - lag;

    if (wh < 0) {
      throw std::range_error("negative wh value");
    } else if (wh < n) {

      out[i] = x[wh];

    } else {
      throw std::range_error("lag is too large");
    }
  }

  return(out);
}



using namespace Rcpp;

//' @title
//' lag_matrix
//'
//' @description
//' lag data and subset the results
//'
//' @inheritParams step_lead_lag
//' @param x to lag (numeric vector)
//' @param lags lead or lag values (numeric vector)
//' @param var_name name for the generated matrix columns (character)
//'
//' @return matrix with lagged values
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericMatrix lag_matrix(const Rcpp::NumericVector& x,
                               const Rcpp::IntegerVector& lags,
                               int n_subset = 1,
                               int n_shift = 0,
                               std::string var_name = "lead_lag") {

  int n = x.size();
  int n_row;

  if(n_subset == 1){
    n_row = (n-n_shift);
  } else {
    n_row = ((n-n_shift-1) / n_subset) + 1;
  }
  int n_col = lags.size();

  Rcpp::CharacterVector nm(n_col);
  Rcpp::NumericMatrix out = Rcpp::NumericMatrix(n_row, n_col);

  for (std::size_t i = 0; i < n_col; i++) {
    out(_, i) = shift_subset(x, lags[i], n_subset, n_shift);
    if(lags[i] < 0) {
      nm[i] = var_name + 'n' + std::to_string(abs(lags[i]));
    } else {
      nm[i] = var_name + std::to_string(lags[i]);
    }
  }

  colnames(out) = nm;

  return(out);
}




struct dl_worker: public Worker {
  // input vector of matrices
  const arma::vec& bv;
  const arma::mat& bl;
  arma::mat& cb;
  int lag_max;
  int n_subset;
  int offset;


  // initialize from Rcpp input and output matrixes (the RMatrix class
  // can be automatically converted to from the Rcpp matrix type)
  dl_worker(const arma::vec& bv,
            const arma::mat& bl,
            arma::mat& cb,
            int lag_max,
            int n_subset,
            int offset)
    : bv(bv), bl(bl), cb(cb), lag_max(lag_max), n_subset(n_subset), offset(offset) {}

  void operator() (std::size_t begin, std::size_t end) {


    for (std::size_t i = begin; i < end; i++) {

      int wh = (i * n_subset) + offset;

      cb.col(i) = bl * bv.subvec(wh, wh + lag_max);

    }
  };
};


//==============================================================================
//' @title
//' distributed_lag_parallel
//'
//' @description
//' This method calculates the basis for a distributed lag in parallel.  It is currently
//' slow.
//'
//' @inheritParams step_lead_lag
//' @param x values to lag (numeric vector)
//' @param bl the basis lags (numeric matrix)
//' @param lag_max maximum number of lags (integer)
//'
//' @return distributed lag basis
//'
//' @noRd
//'
// [[Rcpp::export]]
arma::mat distributed_lag_parallel(const arma::vec& x,
                                   const arma::mat& bl,
                                   int lag_max,
                                   int n_subset = 1,
                                   int n_shift = 0) {

  // result matrix
  int n_row = x.n_elem;
  int n_col = bl.n_rows;

  // added
  //int lag;
  int n_out;
  int start;
  int end;
  int offset;

  if(n_shift >= (n_subset)) {
    throw std::range_error("The absolute value of n_shift should be less than n_subset - 1.");
  }
  if(n_shift < 0) {
    throw std::range_error("n_shift should be positive.");
  }
  if((n_shift + lag_max) > n_row) {
    throw std::range_error("n_shift + lag_max must be less than the length of x");
  }


  //lag   = check_lag(n_row, lag_max, n_shift);
  n_out = get_length(n_row, n_subset);
  Rcpp::NumericVector out(n_out, NA_REAL);

  start = get_start(n_out, lag_max, n_subset);
  end   = get_end(n_row, n_out, lag_max, n_subset);

  if(n_subset != 1) {
    if(n_row % n_subset == 0) {
      offset = n_subset - 1;
    } else {
      offset = abs(n_row - ((n_row / n_subset) * n_subset + 1));
    }
  } else {
    offset = 0;
  }

  int wh = (n_out - start - 1) * n_subset + offset - n_shift + lag_max;

  if (wh < (n_row-n_subset)) {
    start = start - 1;
    offset = offset - n_shift;
  } else {
    offset = offset - n_shift;
  }
  if(offset < 0) {
    n_out = n_out - 1;
    end = end - 1;
    offset = offset + n_subset;
  }


  arma::mat cb(n_col, n_out);
  cb = cb.fill(NA_REAL);

  dl_worker calc_dl(x, bl, cb, lag_max, n_subset, offset);

  RcppParallel::parallelFor(n_out - end, n_out - start, calc_dl);

  return arma::flipud(cb.t());

}


/***R


*/
