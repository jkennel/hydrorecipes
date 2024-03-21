#include "frecipes.h"

//' @title
//' b_spline_list
//'
//' @description
//' Create spline terms
//'
//' @inheritParams splines2::bSpline
//' @param internal_knots locations where parameters can change
//' @param boundary_knots end points of the spline
//' @param complete_basis intercept argument
//'
//' @return List of distributed lags
//'
//' @export
//'
// [[Rcpp::export]]
Rcpp::List b_spline_list(const arma::vec& x,
                         const unsigned int df,
                         const unsigned int degree,
                         const arma::vec& internal_knots,
                         const arma::vec& boundary_knots,
                         const bool complete_basis = false,
                         const bool periodic = false,
                         const unsigned int derivs = 0,
                         const bool integral = false
)
{

  // BSpline object
  splines2::BSpline bs_obj;

  // let splines2 figure out the logic of empty boundary knots
  if (internal_knots.size() > 0) {
    bs_obj = splines2::BSpline(x, internal_knots, degree, boundary_knots);
  } else if (df != 0) {
    bs_obj = splines2::BSpline(x, df, degree, boundary_knots);
  }

  // get B-spline basis functions
  const arma::mat bs_mat = bs_obj.basis(complete_basis);

  size_t n = bs_mat.n_cols;
  Rcpp::List out(n);

  for (size_t i = 0; i < n; ++i) {
    out[i] = bs_mat.col(i);
  }

  return out;
}


//' @title
//' n_spline_list
//'
//' @description
//' Create spline terms
//'
//' @inheritParams splines2::naturalSpline
//' @param internal_knots locations where parameters can change
//' @param boundary_knots end points of the spline
//' @param complete_basis intercept argument
//'
//' @return List of distributed lags
//'
//' @export
//'
// [[Rcpp::export]]
Rcpp::List n_spline_list(const arma::vec& x,
                         const unsigned int df,
                         const unsigned int degree,
                         const arma::vec& internal_knots,
                         const arma::vec& boundary_knots,
                         const bool complete_basis,
                         const bool periodic = false,
                         const unsigned int derivs = 0,
                         const bool integral = false
)
{

  // BSpline object
  splines2::NaturalSpline bs_obj;

  // let splines2 figure out the logic of empty boundary knots
  if (internal_knots.size() > 0) {
    bs_obj = splines2::NaturalSpline(x, internal_knots, boundary_knots);
  } else if (df != 0) {
    bs_obj = splines2::NaturalSpline(x, df, boundary_knots);
  }

  // get natural-spline basis functions
  const arma::mat bs_mat = bs_obj.basis(complete_basis);

  size_t n = bs_mat.n_cols;
  Rcpp::List out(n);

  for (size_t i = 0; i < n; ++i) {
    out[i] = bs_mat.col(i);
  }

  return out;
}



// [[Rcpp::export]]
Rcpp::List b_spline_list2(const arma::vec& x,
                         const unsigned int df,
                         const unsigned int degree,
                         const arma::vec& internal_knots,
                         const arma::vec& boundary_knots,
                         const bool complete_basis = true,
                         const bool periodic = false,
                         const unsigned int derivs = 0,
                         const bool integral = false
)
{

  // BSpline object
  splines2::BSpline bs_obj;

  // let splines2 figure out the logic of empty boundary knots
  if (internal_knots.size() > 0) {
    bs_obj = splines2::BSpline(x, internal_knots, degree, boundary_knots);
  } else if (df != 0) {
    bs_obj = splines2::BSpline(x, df, degree, boundary_knots);
  }

  // get B-spline basis functions
  const arma::mat bs_mat = bs_obj.basis(complete_basis);

  size_t n = bs_mat.n_cols;
  Rcpp::List out(n);

  for (size_t i = 0; i < n; ++i) {
    arma::vec a_vec = bs_mat.col(i);
    Eigen::VectorXd e_vec = Eigen::Map<Eigen::VectorXd>(a_vec.memptr(),
                                                        a_vec.size());
    out[i] = e_vec;
  }

  return out;
}


// [[Rcpp::export]]
std::list<Eigen::VectorXd> b_spline_list3(const arma::vec& x,
                         const unsigned int df,
                         const unsigned int degree,
                         const arma::vec& internal_knots,
                         const arma::vec& boundary_knots,
                         const bool complete_basis = true,
                         const bool periodic = false,
                         const unsigned int derivs = 0,
                         const bool integral = false
)
{

  // BSpline object
  splines2::BSpline bs_obj;

  // let splines2 figure out the logic of empty boundary knots
  if (internal_knots.size() > 0) {
    bs_obj = splines2::BSpline(x, internal_knots, degree, boundary_knots);
  } else if (df != 0) {
    bs_obj = splines2::BSpline(x, df, degree, boundary_knots);
  }

  // get B-spline basis functions
  const arma::mat bs_mat = bs_obj.basis(complete_basis);

  size_t n = bs_mat.n_cols;
  std::list<Eigen::VectorXd> out;
  Eigen::VectorXd tmp_eigen;


  arma::vec tmp_arma(bs_mat.n_cols);
  for (size_t i = 0; i < n; ++i) {
    tmp_arma = bs_mat.col(i);
    tmp_eigen = Eigen::Map<Eigen::VectorXd>(tmp_arma.memptr(), tmp_arma.n_elem);
    out.push_back(tmp_eigen);
  }

  return out;
}

//==============================================================================

//' @title
//' log_lags_arma
//'
//' @description
//' Generate logarithmically spaced lags
//'
//' @param n integer number of lag terms
//' @param max_lag integer the maximum lag
//'
//' @return vector of logarithmically spaced lags
//'
//' @export
//'
// [[Rcpp::export]]
arma::vec log_lags_arma(arma::uword n, arma::uword max_lag) {

  // check inputs
  if (n <= 0) {
    Rcpp::stop("log_lags_eigen: n must be greater than 0");
  }

  if (max_lag < 0) {
    Rcpp::stop("log_lags_eigen: max_time_lag must be non-negative");
  }

  if (n > (max_lag + 1L)) {
    Rcpp::warning("The number of lags is greater than the maximum time lag");
    return(arma::linspace(0.0, (double)max_lag, max_lag + 1));
  }
  // Lags begin at zero
  arma::vec lags = arma::exp(arma::linspace(0.0, std::log((double)max_lag + 1.0), n))-1;

  // lags cannot be spaced closer than one sample
  for (size_t i = 0; i < n; ++i) {
    if (lags(i) < i) {
      lags(i) = i;
    }
  }

  return(round(lags));
}

/*** R

n <- 2e6
m <- sort(rnorm(n))
bk <- range(m)
knots <- quantile(bk, probs = seq(0.05, 0.95, 0.3))
bench::mark(
  tmp <- frecipes:::b_spline_list(m, 0L, 3L, knots, bk),
  tmp <- frecipes:::b_spline_list2(m, 0L, 3L, knots, bk),
  tmp <- frecipes:::b_spline_list3(m, 0L, 3L, knots, bk),
  check = FALSE,
  min_iterations = 5
)

microbenchmark::microbenchmark(
  # tmp <- rcpp_bSpline_fit(m, 0L, 3L, knots, bk),
  tmp3 <- b_spline(m, c(bk[1], knots, bk[2]), 3),
  times = 2
)

bench::mark(
  hydrorecipes:::log_lags_eigen(1000, max_lag = 1e7),
  log_lags_arma(1000, max_lag = 1e7),
  check = TRUE
)
*/
