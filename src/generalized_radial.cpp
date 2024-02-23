#include "frecipes.h"


//==============================================================================

// [[Rcpp::export]]
std::vector<double> impulse_function(std::vector<double> u)
{

  int n = u.size();

  // calculate the pulse
  for (unsigned int i = n - 1; i > 0; --i) {
    u[i] -= u[i-1];
  }

  return u;
}

// [[Rcpp::export]]
Rcpp::NumericVector impulse_function_rcpp(Rcpp::NumericVector u)
{

  int n = u.size();

  // calculate the pulse
  for (size_t i = n - 1; i > 0; --i) {
    u[i] -= u[i-1];
  }

  return u;
}

// [[Rcpp::export]]
Eigen::VectorXd impulse_function_eigen(Eigen::VectorXd u)
{

  int n = u.size();

  // calculate the pulse
  for (size_t i = n - 1; i > 0; --i) {
    u[i] -= u[i-1];
  }

  return u;
}

//==============================================================================

// [[Rcpp::export]]
double exp_int(double u) {

    if (u == 0){
      u = R_PosInf;
    } else if (u > 700.0){
      u = 0;
    } else {
      u = -std::expint(-u);
    }

  return(u);
}

// [[Rcpp::export]]
int binary_search(Eigen::VectorXd x, Eigen::VectorXd y)
{
    int r = x.size() - 1;
    int l = 0;

    if (x[l] >= y[l]) {
      return(l);
    }

    if (x[r] < y[r]) {
      return(r);
    }

    while (l <= r) {
      int m = l + (r - l) / 2;

      // If x greater, ignore right half
      if (x[m] >= y[m]) {
        if (x[m-1] < y[m-1]) {
          return(m);
        } else {
          r = m - 1;
          // If x is smaller, ignore left half
        }} else {
          l = m + 1;
        }

    }

    // unsuccessful
    return -1;
}

// [[Rcpp::export]]
Eigen::VectorXd std_to_eigen(std::vector<double> u)
{
  Eigen::VectorXd out = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(u.data(), u.size());
  return(out);
}

// [[Rcpp::export]]
std::vector<double> eigen_to_std(Eigen::VectorXd u)
{
  std::vector<double> out(u.data(), u.data() + u.size());
  return(out);
}

// [[Rcpp::export]]
Rcpp::NumericVector std_to_rcpp(std::vector<double> u)
{
  return(Rcpp::NumericVector(u.begin(), u.end()));
}

// [[Rcpp::export]]
std::vector<double> rcpp_to_std(Rcpp::NumericVector u)
{
  return(Rcpp::as<std::vector<double> >(u));
}


//==============================================================================

// [[Rcpp::export]]
double calculate_distance(const double x_well,
                          const double y_well,
                          const double x_loc,
                          const double y_loc
)
{

  double dx = x_well - x_loc;
  double dy = y_well - y_loc;
  return(std::sqrt(dx * dx + dy * dy));

}



//==============================================================================

// [[Rcpp::export]]
double well_function_coefficient(const double flow_rate,
                                 const double transmissivity)
{
  return(flow_rate / (4.0 * M_PI * transmissivity));
}


// [[Rcpp::export]]
std::vector<double> well_function_coefficient_vec(std::vector<double> flow_rate,
                                                  const double transmissivity)
{
  for (auto &out : flow_rate)
    out = well_function_coefficient(out, transmissivity);

  return flow_rate;
}

// [[Rcpp::export]]
Rcpp::NumericVector well_function_coefficient_rcpp(Rcpp::NumericVector flow_rate,
                                                       const double transmissivity)
{
  return flow_rate / (4.0 * M_PI * transmissivity);
}




// [[Rcpp::export]]
double theis_u(double radius,
               double storativity,
               double transmissivity,
               double time) {

  return ((radius * radius * storativity) / (4.0 * transmissivity * time));

}

// [[Rcpp::export]]
std::vector<double> theis_u_time_vec(double radius,
                                     double storativity,
                                     double transmissivity,
                                     std::vector<double> time) {


  for (auto &out : time)
    out = theis_u(radius, storativity, transmissivity, out);

  return time;
}

// [[Rcpp::export]]
Rcpp::NumericVector theis_u_time_rcpp(const double radius,
                                      const double storativity,
                                      const double transmissivity,
                                      const Rcpp::NumericVector time) {

  return ((radius * radius * storativity) / (4.0 * transmissivity * time));
}


// [[Rcpp::export]]
double grf_coefficient(const double radius,
                       const double hydraulic_conductivity,
                       const double thickness,
                       const double flow_dimension)
{

  double v = 1.0 - flow_dimension / 2.0;

  return ((pow(radius, (2.0 * v))) /
          (4.0 * pow(M_PI, (1.0 - v)) * hydraulic_conductivity *
            pow(thickness, (3.0 - flow_dimension))));
}






// Time is dealt with later
// [[Rcpp::export]]
double grf_u(const double radius,
             const double specific_storage,
             const double hydraulic_conductivity)
{

  return ((radius * radius * specific_storage) /
          (4.0 * hydraulic_conductivity));
  ;
}


//==============================================================================
//' @title
//' grf_time
 //'
 //' @description
 //' Parallel convolution of GRF well function and flow rates in the time domain.
 //' Time series needs to be regularily spaced and so are the flow rates.  Some
 //' performance gains can be achieved if the number of flow rate does not change
 //' for each time.
 //'
 //' @param radius distance to monitoring interval
 //' @param specific_storage aquifer storativity
 //' @param hydraulic_conductivity aquifer hydraulic conductivity
 //' @param thickness aquifer thickness
 //' @param time prediction times
 //' @param flow_rate well flow rates
 //' @param flow_time_interval time between flow rate measurements in samples
 //' @param flow_dimension flow dimension
 //'
 //' @return theis solution for multiple pumping scenario
 //'
 //'
 //' @export
 //'
// [[Rcpp::export]]
Rcpp::List grf_time(const double radius,
                    const double specific_storage,
                    const double hydraulic_conductivity,
                    const double thickness,
                    const Rcpp::NumericVector time,
                    const Rcpp::NumericVector flow_rate,
                    const double flow_dimension)
{

  size_t n_flow_rate = flow_rate.size();
  size_t n_time = time.size();

  // check that the number of times and flow rates are equal
  if (n_flow_rate != n_time)
  {
    Rcpp::stop("The number of times and flow_rate should be the same");
  }

  const double v = (flow_dimension / 2.0) - 1.0;

  // calculate the constant part
  double u_const = grf_u(radius, specific_storage, hydraulic_conductivity);
  double coef_const = grf_coefficient(radius,
                                      hydraulic_conductivity,
                                      thickness,
                                      flow_dimension);

  Rcpp::NumericVector coef = coef_const * flow_rate;
  Rcpp::NumericVector u = u_const / time;

  u = specialfunctions::gamma_inc_rcpp(u, v);
  u = impulse_function_rcpp(u);

  Eigen::VectorXd u_eig(Rcpp::as<Eigen::VectorXd>(u));
  Eigen::VectorXd coef_eig(Rcpp::as<Eigen::VectorXd>(coef));;

  return Rcpp::List::create(
    Rcpp::Named("generalized_radial") = convolve_filter(u_eig, coef_eig, false, true)
  );

}


//==============================================================================
//' @title
 //' grf_grid
 //'
 //' @description
 //' Parallel convolution of GRF well function and flow rates in the time domain.
 //' Time series needs to be regularily spaced and so are the flow rates.  Some
 //' performance gains can be achieved if the number of flow rate does not change
 //' for each time.
 //'
 //' @param radius distance to monitoring interval
 //' @param specific_storage aquifer storativity
 //' @param hydraulic_conductivity aquifer hydraulic conductivity
 //' @param thickness aquifer thickness
 //' @param time prediction times
 //' @param flow_rate well flow rates
 //' @param flow_time_interval time between flow rate measurements in samples
 //' @param flow_dimension flow dimension
 //'
 //' @return theis solution for multiple pumping scenario
 //'
 //'
 //' @export
 //'
 // [[Rcpp::export]]
 Eigen::MatrixXd grf_grid(const Eigen::MatrixXd &grid,
                          const Eigen::MatrixXd &well_locations,
                          const Eigen::MatrixXd &flow_rate,
                          const Eigen::VectorXd &time,
                          const double specific_storage,
                          const double hydraulic_conductivity,
                          const double thickness,
                          const double flow_dimension)
 {

   size_t n_grid = grid.rows();
   size_t n_well = well_locations.rows();

   size_t n_flow_rate = flow_rate.rows();
   size_t n_flow_rate_well = flow_rate.cols();
   size_t n_time = time.size();

   // check that the number of times and flow rates are equal
   if (n_flow_rate != n_time)
   {
     Rcpp::stop("The number of times and flow_rate should be the same");
   }
   if (n_well != n_flow_rate_well)
   {
     Rcpp::stop("Dimensions of flow_rate and well_locations should be consistent");
   }
   const double v = (flow_dimension / 2.0) - 1.0;

   // calculate distances

   Eigen::MatrixXd output = Eigen::MatrixXd::Zero(n_grid, n_time);

   // RcppThread::parallelFor (0, n_grid, [&] (size_t i) {

   Eigen::VectorXd coef(n_time);
   Eigen::VectorXd wf(n_time);
   // std::vector<double> impulse(n_time);
   std::vector<double> u(n_time);
   double distance;
   double u_const;
   double coef_const;

   for (size_t i = 0; i < n_grid; ++i) {
     for (size_t j = 0; j < n_well; ++j) {

       distance = calculate_distance(
         grid(i, 0),
         grid(i, 1),
         well_locations(j, 0),
         well_locations(j, 1));
       u_const = grf_u(distance, specific_storage, hydraulic_conductivity);
       coef_const = grf_coefficient(distance,
                                    hydraulic_conductivity,
                                    thickness,
                                    flow_dimension);
       coef = coef_const * flow_rate.col(j);
       wf = u_const / time.array();
       VectorXd::Map(&u[0], n_time) = wf;

       u = specialfunctions::gamma_inc_vec(u, v);
       u = impulse_function(u);

       wf = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(u.data(), n_time);

       output.row(i) += convolve_filter(wf, coef, false, true);

     }
   }

   return (output);
 }



//==============================================================================


// [[Rcpp::export]]
double hantush_epsilon(const double radius,
                       const double leakage)
{
  return((radius * radius) / (4.0 * leakage * leakage));
}


//==============================================================================
//' @title
//' hantush_well_single
//'
//' @description
//' Result of the hantush well function
//'
//' J.H.A. Prodanoff; W.J. Mansur; F.C.B. Mascarenhas (2006). Numerical
//' evaluation of Theis and Hantush-Jacob well functions. , 318(1-4),
//' 0–183. doi:10.1016/j.jhydrol.2005.05.026 eq: 10, 11, 12
//'
//'
//' @param u value of the Theis u
//' @param b the leakance
//' @param n_terms the number of terms used in the hantush approximation
//'
//'
//' @return hantush well function
//'
//'
//' @export
//'
// [[Rcpp::export]]
double hantush_well(double u, double b, double precision){

  double b_div_u = b / u;
  double out = 0.0;
  double en;
  double to_add;
  unsigned int n_terms = 30;

  //eq 10
  if (b_div_u >= u){
    en = exp_int(b_div_u);

    for (unsigned int i = 0; i < n_terms; i++) {
      to_add = en * (pow(-u, i) / std::tgamma(i+1));
      out += to_add;
      if (std::abs(to_add) < precision){
        break;
      }
      en = (1.0 / ((double)i + 1.0)) * (exp(-b_div_u) - b_div_u * en);

    }
    out = 2.0 * std::cyl_bessel_k(0, 2.0 * sqrt(b)) - out;

  } else { //eq 12

    en = exp_int(u);

    for (unsigned int i = 0; i < n_terms; i++) {
      to_add = en * (pow(-b_div_u, i) / std::tgamma(i + 1));
      out += to_add;
      if (std::abs(to_add) < precision){
        break;
      }
      en = (1.0 / ((double)i + 1.0)) * (exp(-u) - u * en);

    }

  }

  if (out < 0.0) {
    out = 0.0;
  }

  return(out);
}

// [[Rcpp::export]]
std::vector<double> hantush_well_vec(std::vector<double> u, double b, int n_terms){

  for (auto& out : u)
    out = hantush_well(out, b, n_terms);

  return(u);
}

// [[Rcpp::export]]
Rcpp::NumericVector hantush_well_rcpp(Rcpp::NumericVector u, double b, double precision){

  unsigned int n = u.size();
  Rcpp::NumericVector out(n);

  for (unsigned int i = 0; i < n; ++i){
    out(i) = hantush_well(u(i), b, precision);
  }

  return(out);
}



//==============================================================================
//' @title
//' hantush_jacob
//'
//' @description
//' Convolution of hantush well function and flow rates in the time domain.
//' Time series needs to be regularly spaced.
//'
//' @param radius distance to monitoring interval
//' @param storativity aquifer storativity
//' @param transmissivity aquifer transmissivity
//' @param leakage hantush leakage
//' @param time prediction times
//' @param flow_rate well flow rates
//' @param flow_time_interval time between flow rate measurements in samples
//' @param n_terms number of terms to use in Hantush solution.  More is more precise but slower.
//'
//' @return hantush jacob solution for multiple pumping scenario
//'
//'
//' @export
//'
// [[Rcpp::export]]
Rcpp::List hantush_jacob(
    const Rcpp::NumericVector time,
    const Rcpp::NumericVector flow_rate,
    const double radius,
    const double storativity,
    const double transmissivity,
    const double leakage,
    const double precision
) {
  const int n = time.size();
  const double b = hantush_epsilon(radius, leakage);

  Rcpp::NumericVector coef = well_function_coefficient_rcpp(flow_rate, transmissivity);

  Rcpp::NumericVector u = theis_u_time_rcpp(radius, storativity, transmissivity, time);

  u = hantush_well_rcpp(u, b, precision);

  u = impulse_function_rcpp(u);

  Eigen::VectorXd u_eig(Rcpp::as<Eigen::VectorXd>(u)); //Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(u.data(), u.size());
  Eigen::VectorXd coef_eig(Rcpp::as<Eigen::VectorXd>(coef));// = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(u.data(), u.size());

  return Rcpp::List::create(
    Rcpp::Named("hantush_jacob") = convolve_filter(u_eig, coef_eig, false, true)
  );

}




/*** R
y <- rev(sort(abs(rnorm(1000000))))
x <- 1/y
bench::mark(
frecipes:::binary_search(x, y)
)

x <- rnorm(1e6)
y <- 1.2
bench::mark(
frecipes:::well_function_coefficient_vec_rcpp(x,y),
frecipes:::well_function_coefficient_vec(x,y),
check = FALSE
)

x <- abs(rnorm(1000000))
bench::mark(frecipes:::ei_bh_vec(x),
            frecipes:::ei_sp_vec(x),
            frecipes:::ei_vec(x),
            frecipes:::ei_eigen(x),
            iterations = 1)

n <- 1e5
storativity = 1e-5
radius = 50
transmissivity = 1e-3
leakage <- 100
times <- as.numeric(1:n)
flow_rates <- rep(1, n)#abs(rnorm(n))
n_terms <- 12L

bench::mark(

  tmp <- frecipes:::hantush_jacob(
    times,
    flow_rates,
    radius,
    storativity,
    transmissivity,
    leakage,
    n_terms)[[1]],

  tmp2 <- aquifer:::hantush_convolve(radius,
                                     storativity,
                                     transmissivity,
                                     leakage,
                                     times,
                                     flow_rates,
                                     n_terms),
  check = FALSE
)

# x <- rnorm(100000)
# bench::mark(frecipes:::convolve_overlap_add(x,x),
#             frecipes:::convolve_filter(x,x,TRUE, TRUE), check = FALSE)

plot(tmp2[1:1000], col = 'red', type = 'l', log = 'xy')
points(tmp[1:1000], type = 'l', log = 'xy')

x <- sort(abs(rnorm(10000)))
bench::mark(
  # frecipes:::hantush_well_vec(x, 0.01, 10),
  # frecipes:::hantush_well_e(x, 0.01, 10),
  frecipes:::hantush_well_rcpp(x, 0.01, 10),
  aquifer:::hantush_well_parallel(x, 0.01, 10),
  check = FALSE
)


bench::mark(
  frecipes:::hantush_jacob(x, rep(0.01, length(x)),
                           radius = 10, 1e-6, 1e-3, 1, 20),
  sapply(1:100, function(x) aquifer:::hantush_well_single(0.1, 0.01, 10)),
  check = FALSE
)

x <- abs(rnorm(10000))
bench::mark(
  frecipes:::bessel_eigen(x),
  frecipes:::bessel_rcpp(x),
  besselK(x, 1),
  Bessel::BesselK(x, 1)
)

library(expint)
x <- abs(rnorm(10000))
bench::mark(frecipes:::gis(x, 0.0),
            frecipes:::gis2(x),
            check = FALSE)

n <- 1000000
time <- seq(1, n, 1)
flow_rate <- rep(0.001, n)
thickness <- 10
radius <- 5
specific_storage_1 <- 1e-5
hydraulic_conductivity_1 <- 1e-5
specific_storage_2 <- 5e-6
hydraulic_conductivity_2 <- 5e-6
diffusivity_1 <- hydraulic_conductivity_1 / specific_storage_1
diffusivity_2 <- hydraulic_conductivity_2 / specific_storage_2
bench::mark(

grf_1 <- frecipes:::grf_time(radius,
                             specific_storage_1,
                             hydraulic_conductivity_1,
                             thickness,
                             time,
                             flow_rate,
                             flow_dimension = 2)
)
plot(grf_1[[1]], type = 'l', log = 'xy')
points(grf_1[[1]], type = 'l', col = 'red')
points(grf_1[[1]], type = 'l', col = 'blue')
specific_storage <- 0.5e-5
hydraulic_conductivity <- 0.5e-5
grf_2 <- grf_time(radius,
                  specific_storage_2,
                  hydraulic_conductivity_2,
                  thickness,
                  time,
                  flow_rate,
                  flow_dimension = 2)

dat <- data.table(time, grf_1, grf_2)
fit_1 <- lm(grf_1~log(time), tail(dat, 4000))
fit_2 <- lm(grf_2~log(time), tail(dat, 4000))
summary(fit_1)
summary(fit_2)
plot(grf_2, log = 'x', type = 'l', xlab = "Elapsed time", ylab = "drawdown")
abline(h = 0, col = 'grey')
points(grf_1, type = 'l', col = 'red')

dat[, pred_1 := predict(fit_1, dat)]
dat[, pred_2 := predict(fit_2, dat)]

points(pred_1~time, dat, type = 'l', col = 'red', lty = 2)
points(pred_2~time, dat, type = 'l', lty = 2)

*/
