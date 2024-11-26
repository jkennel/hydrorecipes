#include "hydrorecipes.h"


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
double std_expint(double u) {

    if (u == 0){
      u = R_PosInf;
    // } else if (u > 708.0){
    //   u = 708.0;
    } else {
      u = -std::expint(-u);
    }

  return(u);
}

// [[Rcpp::export]]
Eigen::VectorXd std_expint_vec(std::vector<double> u) {

  for (auto &out : u)
    out = std_expint(out);


  return(std_to_eigen(u));
}


// [[Rcpp::export]]
double std_tgamma(double u, double a) {

  double ret;

  if (std::isfinite(u)) {
    ret = std::tgamma(u);
  }
  else {
    ret = 0.0;
  }

  return(ret);
}

// [[Rcpp::export]]
double bh_gamma_p_inv(double a, double p) {

  return(boost::math::gamma_p_inv(a, p));

}


//
//
// // [[Rcpp::export]]
// double bh_gamma_neg(double u, double a) {
//
//   double ret;
//
//   if (std::isfinite(u)) {
//     ret = ((1.0 - boost::math::gamma_p(a + 1.0, u)) * std::tgamma(a + 1.0) -
//       std::pow(u, a) * std::exp(-u)) / a;
//   }
//   else {
//     ret = 0.0;
//   }
//
//   return(ret);
//
// }
//
// // [[Rcpp::export]]
// double gamma_der(double u, double a) {
//   double ret;
//
//   if (std::isfinite(u)) {
//     ret = std::pow(u, a - 1) / std::exp(u);
//   }
//   else {
//     ret = 0.0;
//   }
//   return(ret);
// }
//
//
// // [[Rcpp::export]]
// double gamma_inc(double u, double a) {
//   double ret;
//
//   if(a==0){
//     ret = std_exp_int(u);
//   } else if(a>0){
//     ret = std_tgamma(u, a);
//   } else {
//     ret = std_gamma_neg(u, a);
//   }
//
//
//   return(ret);
// }







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
double theis_u(const double radius,
               const double storativity,
               const double transmissivity,
               const double time) {

  return ((radius * radius * storativity) / (4.0 * transmissivity * time));

}


// [[Rcpp::export]]
std::vector<double> theis_u_time_vec(const double radius,
                                     const double storativity,
                                     const double transmissivity,
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
double theis_aniso_coefficient(const double transmissivity_x,
                               const double transmissivity_y) {

  return(1.0 / (4.0 * M_PI * std::sqrt(transmissivity_x * transmissivity_y)));

}

// [[Rcpp::export]]
double theis_aniso_u(const double x,
                     const double y,
                     const double storativity,
                     const double transmissivity_x,
                     const double transmissivity_y)
{

  return ( storativity / (4.0) *
           (transmissivity_x * y * y + transmissivity_y * x * x) /
           (transmissivity_x * transmissivity_y));
}

// [[Rcpp::export]]
Eigen::VectorXd theis_aniso_u_grid(Eigen::VectorXd x,
                                   Eigen::VectorXd y,
                                   const double storativity,
                                   const double transmissivity_x,
                                   const double transmissivity_y)
{

  return (storativity / (4.0) *
          (transmissivity_x * y.array().square() + transmissivity_y * x.array().square()) /
            (transmissivity_x * transmissivity_y));
}


// // [[Rcpp::export]]
// std::vector<double> theis_aniso_u_time_vec(const double x,
//                                            const double y,
//                                            const double storativity,
//                                            const double transmissivity_x,
//                                            const double transmissivity_y,
//                                            std::vector<double> time) {
//
//
//   for (auto &out : time)
//     out = theis_aniso_u(x, y, storativity, transmissivity_x, transmissivity_y, out);
//
//   return time;
// }


//==============================================================================
//' @title
//' theis_aniso_time
//'
//' @description
//' Convolution of GRF well function and flow rates in the time domain.
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
//' @noRd
// [[Rcpp::export]]
Rcpp::List theis_aniso_time(const double distance_x,
                            const double distance_y,
                            const double storativity,
                            const double transmissivity_x,
                            const double transmissivity_y,
                            const double thickness,
                            Eigen::VectorXd time,
                            Eigen::VectorXd flow_rate)
{

  size_t n_flow_rate = flow_rate.size();
  size_t n_time = time.size();

  // check that the number of times and flow rates are equal
  if (n_flow_rate != n_time)
  {
    Rcpp::stop("The number of times and flow_rate should be the same");
  }

  const double a = -1.0;

  // calculate the constant part
  const double u_const = theis_aniso_u(distance_x, distance_y,
                                       storativity,
                                       transmissivity_x,
                                       transmissivity_y);

  const double coef_const = theis_aniso_coefficient(transmissivity_x,
                                                    transmissivity_y);

   Eigen::VectorXd coef = coef_const * flow_rate.array();
   Eigen::VectorXd u = u_const / time.array();

   u = impulse_function_eigen(std_expint_vec(eigen_to_std(u)));

   return Rcpp::List::create(
     Rcpp::Named("generalized_radial") = convolve_filter(u, coef, false, true)
   );

 }



// [[Rcpp::export]]
Rcpp::List grid_pumping_regimes(Eigen::VectorXd distance_x,
                                Eigen::VectorXd distance_y,
                                Eigen::VectorXd output_times,
                                Eigen::VectorXd start_times,
                                Eigen::VectorXd flow_rates,
                                Eigen::VectorXd well_x,
                                Eigen::VectorXd well_y,
                                double storativity,
                                double transmissivity_x,
                                double transmissivity_y,
                                double thickness)
{

  unsigned int n_out_times = output_times.size();
  unsigned int n_in_times = start_times.size();
  unsigned int n_well = well_x.size();
  unsigned int n_obs = distance_x.size();

  Eigen::VectorXd dx = Eigen::ArrayXd::Zero(n_obs);
  Eigen::VectorXd dy = dx;

  Eigen::VectorXd u_const(n_obs);

  const double coef_const = theis_aniso_coefficient(transmissivity_x, transmissivity_y);

  Rcpp::List out(n_out_times);
  for (unsigned int i = 0; i < n_out_times; ++i) {
    out[i] = Eigen::VectorXd::Zero(n_obs);
  }

  for (unsigned int j = 0; j < n_in_times; ++j) {
    bool new_well = (j == 0) || (well_x(j - 1) != well_x(j) || well_y(j - 1) != well_y(j));
    double flow_change = new_well ? flow_rates(j) : (flow_rates(j) - flow_rates(j - 1));

    if (new_well) {
      u_const = theis_aniso_u_grid(distance_x.array() - well_x(j),
                                   distance_y.array() - well_y(j),
                                   storativity,
                                   transmissivity_x,
                                   transmissivity_y).array();
    }

    const double coef = coef_const * flow_change;

    for (unsigned int i = 0; i < n_out_times; ++i) {
      double elapsed_time = output_times(i) - start_times(j);
      if (elapsed_time > 0) {
        out[i] = Rcpp::as<Eigen::VectorXd>(out[i]).array() +
          coef * std_expint_vec(eigen_to_std(u_const.array() / elapsed_time)).array();
      }
    }
  }

  return out;
}





// // well_rates
// // 0) x
// // 1) y
// // 2) time
// // 3) flow_rate
// // [[Rcpp::export]]
// Rcpp::List theis_aniso_grid_by_well(Eigen::VectorXd monitor_x,
//                                          Eigen::VectorXd monitor_y,
//                                          Eigen::MatrixXd well_rates,
//                                          const double storativity,
//                                          const double transmissivity_x,
//                                          const double transmissivity_y,
//                                          const double thickness,
//                                          Eigen::VectorXd output_times)
// {
//
//   const unsigned int n_well = well_rates.rows();
//   const unsigned int n_monitor = monitor_x.size();
//   const unsigned int n_times = output_times.size();
//
//   Eigen::VectorXd well_contrib = Eigen::VectorXd::Zero(n_monitor);
//   Rcpp::List out(n_times);
//   Eigen::VectorXd distance_x(n_monitor);
//   Eigen::VectorXd distance_y(n_monitor);
//
//   double flow_rate;
//   double output_time;
//   double elapsed_time;
//
//   for (unsigned int j = 0; j < n_times; ++j) {
//
//     output_time = output_times(j);
//
//     for (unsigned int i = 0; i < n_well; ++i) {
//
//       elapsed_time = output_time - well_rates(i, 2);
//
//       if (elapsed_time > 0) {
//
//         distance_x = monitor_x.array() - (double)well_rates(i, 0);
//         distance_y = monitor_y.array() - (double)well_rates(i, 1);
//         flow_rate = 1.0;
//
//         well_contrib += theis_aniso_grid(distance_x,
//                 distance_y,
//                 storativity,
//                 transmissivity_x,
//                 transmissivity_y,
//                 thickness,
//                 elapsed_time,
//                 flow_rate);
//       }
//
//       out[j] = well_contrib;
//
//     }
//
//
//   }
//
//   return(out);
// }



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
//' Convolution of GRF well function and flow rates in the time domain.
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
//' @noRd
// [[Rcpp::export]]
Rcpp::List grf_time(const double radius,
                    const double specific_storage,
                    const double hydraulic_conductivity,
                    const double thickness,
                    Eigen::VectorXd time,
                    Eigen::VectorXd flow_rate,
                    const double flow_dimension)
{

  size_t n_flow_rate = flow_rate.size();
  size_t n_time = time.size();

  // check that the number of times and flow rates are equal
  if (n_flow_rate != n_time)
  {
    Rcpp::stop("The number of times and flow_rate should be the same");
  }

  double a = (flow_dimension / 2.0) - 1.0;

  // calculate the constant part
  double u_const = grf_u(radius, specific_storage, hydraulic_conductivity);
  double coef_const = grf_coefficient(radius,
                                      hydraulic_conductivity,
                                      thickness,
                                      flow_dimension);

  Eigen::VectorXd coef = coef_const * flow_rate.array();
  Eigen::VectorXd u = u_const / time.array();

  u = gamma_inc(u.array(), a);
  u = impulse_function_eigen(u);

  return Rcpp::List::create(
    Rcpp::Named("generalized_radial") = convolve_filter(u, coef, false, true)
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
//' @noRd
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
  // std::vector<double> u(n_time);
  Eigen::VectorXd u(n_time);

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

      u = gamma_inc(u.array(), v);
      // VectorXd::Map(&u[0], n_time) = wf;

      // u = specialfunctions::gamma_inc_vec(u, v);
      u = impulse_function_eigen(u);

      // wf = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(u.data(), n_time);

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
//' hantush_well
//'
//' @description
//' Result of the hantush well function
//'
//' Prodanoff, J.H.A., Mansur, W.J. and Mascarenhas, F.C.B., 2006. Numerical
//'   evaluation of Theis and Hantush-Jacob well functions. Journal of
//'   hydrology, 318(1-4), pp.173-183.
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
//' @noRd
//'
// [[Rcpp::export]]
double hantush_well(double u, double b, double precision){

  double b_div_u = b / u;
  double out = 0.0;
  double en;
  double to_add;
  unsigned int n_terms = 64;

  //eq 10
  if (b_div_u >= u){
    // Rcpp::Rcout << "here: " << b_div_u << std::endl;
    en = std_expint(b_div_u);

    for (unsigned int i = 0; i < n_terms; ++i) {
      to_add = en * (pow(-u, i) / std::tgamma(i + 1));
      out += to_add;
      if (std::fabs(to_add) < precision){
        break;
      }
      en = (1.0 / ((double)i + 1.0)) * (exp(-b_div_u) - b_div_u * en);

    }
    out = 2.0 * std::cyl_bessel_k(0, 2.0 * sqrt(b)) - out;

  } else { //eq 12
    // Rcpp::Rcout << "there: " << b_div_u << std::endl;

    en = std_expint(u);

    for (unsigned int i = 0; i < n_terms; ++i) {
      to_add = en * (pow(-b_div_u, i) / std::tgamma(i + 1)); // tgamma(i+1) = factorial(i)
      out += to_add;
      if (std::fabs(to_add) < precision){
        break;
      }
      en = (1.0 / ((double)i + 1.0)) * (exp(-u) - u * en);
    }

  }
  // Rcpp::Rcout << "out: " << out << std::endl;

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
//' @noRd
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

// [[Rcpp::export]]
Eigen::VectorXd ig(Eigen::ArrayXd a, Eigen::ArrayXd u) {
  return(Eigen::igammac(a, u));
}



// (https://en.wikipedia.org/wiki/Exponential_distribution)
// lamda is a rate parameter
// amplitude is a scaling parameter
// [[Rcpp::export]]
double exp_2_parameter(double t,
                       const double amplitude,
                       const double lamda) {

  return (amplitude / lamda) * std::exp(-t / lamda);

}

// [[Rcpp::export]]
std::vector<double> exp_2_old(std::vector<double> t, double amplitude, double lamda){


  for (auto& out : t)
    out = exp_2_parameter(out, amplitude, lamda);

  return(t);
}

// [[Rcpp::export]]
std::vector<double> exp_2(std::vector<double> t, const double amplitude, const double lamda){

  const double scale = amplitude / lamda;

  double exp_val;

  for (auto& out : t) {
    exp_val = out / -lamda;
    if (exp_val < -800) {
      out = 0;
    } else {
      out = scale * std::exp(exp_val);
    }
  }
  return(t);
}

// [[Rcpp::export]]
std::vector<double> exp_2_test(const double amplitude,
                               const double lamda){

  const double scale = amplitude / lamda;
  int max_t = 750 * std::ceil(lamda);

  std::vector<double> out(max_t);
  for (int i = 0; i < max_t; ++i) {
    out[i] = scale * std::exp((double)i / -lamda);
  }
  return(out);
}


// [[Rcpp::export]]
Eigen::ArrayXd exp_2_eigen(const double amplitude,
                           const double lamda){

  const int max_t = 720 * std::ceil(lamda);

  Eigen::ArrayXd out = Eigen::ArrayXd::LinSpaced(max_t + 1, 0.0, (double)max_t);

  return(amplitude / lamda * (out / -lamda).exp());

}



// for specifying response functions
// (https://en.wikipedia.org/wiki/Gamma_distribution)
// k is a shape parameter
// theta is a scaling parameter
// amplitude is a scaling parameter
// [[Rcpp::export]]
double gamma_3_parameter(double t,
                         const double amplitude,
                         const double k,
                         const double theta) {

  return (amplitude * std::pow(t,  (k - 1.0)) * std::exp(-t / theta) /
          (std::pow(theta, k) * std::tgamma(k)));

}

// [[Rcpp::export]]
std::vector<double> gamma_3_old(std::vector<double> t,
                                const double amplitude,
                                const double k,
                                const double theta){


  for (auto& out : t)
    out = gamma_3_parameter(out, amplitude, k, theta);

  return(t);
}


// [[Rcpp::export]]
std::vector<double> gamma_3_old2(std::vector<double> t,
                                 const double amplitude,
                                 const double k,
                                 const double theta){

  const double denom = (std::pow(theta, k) * std::tgamma(k));

  for (auto& out : t)
    out = amplitude * std::pow(out,  (k - 1.0)) * std::exp(out / -theta) / denom;

  return(t);
}

// [[Rcpp::export]]
std::vector<double> gamma_3(std::vector<double> t,
                            const double amplitude,
                            const double k,
                            const double theta){


  // return exponetial if k = 1
  if (k == 1) {
    return(exp_2(t, amplitude, theta));
  }

  const double denom = (std::pow(theta, k) * std::tgamma(k));
  double exp_val;

  for (auto& out : t) {
    exp_val = out / -theta;
    if (exp_val < -800) {
      out = 0;
    } else {
      out = amplitude * std::pow(out, (k - 1.0)) * std::exp(exp_val) / denom;
    }
  }


  return(t);
}

// [[Rcpp::export]]
std::vector<double> gamma_3_test(const double amplitude,
                                 const double k,
                                 const double theta){

  const double denom = (std::pow(theta, k) * std::tgamma(k));
  double exp_val;
  int max_t = 720 * std::ceil(theta);

  std::vector<double> out(max_t + 1);
  for (int i = 0; i < max_t + 1; ++i) {
    out[i] = amplitude * std::pow(i, (k - 1.0)) * std::exp(((double)i / -theta)) / denom;
  }

  return(out);
}

// [[Rcpp::export]]
Eigen::ArrayXd gamma_3_eigen(const double amplitude,
                             const double k,
                             const double theta){

  const double denom = (std::pow(theta, k) * std::tgamma(k));
  const int max_t = 720 * std::ceil(theta);

  Eigen::ArrayXd out = Eigen::ArrayXd::LinSpaced(max_t + 1, 0.0, (double)max_t);

  return(amplitude * out.pow(k - 1.0) * (out / -theta).exp() / denom);

  // for (int i = 1; i < max_t; ++i) {
  //   out[i] = amplitude * std::pow(i, (k - 1.0)) * std::exp(((double)i / -theta)) / denom;
  // }
  //
  // return(out);
}

// // [[Rcpp::export]]
// std::vector<double> test(std::vector<double> t) {
//
//   std::gamma_distribution<double> d(1.0, 2.0);
//
//   for (auto& out : t)
//     out = d(out);
//
//   return(t);
//
// }

// generate all u (each grid point each well)
// determine range
// generate vector to cover u values (intelligent spacing)
// Eigen::VectorXd interpolate_grid(Eigen::VectorXd x_in,
//                             Eigen::VectorXd y_in,
//                             Eigen::VectorXd x_out) {
//
//   unsigned int n = x_in.size();
//   Eigen::VectorXd y_out(n);
//   return(y_out);
//
// }

/*** R

library(bench)
library(hydrorecipes)

A <- 1.4
a <- 10000
n <- 1.0
t <- as.numeric(1:1e6)


bench::mark(
hydrorecipes:::bh_gamma_p_inv(n, 0.9999) * a
)

tmp <- hydrorecipes:::gamma_3(t, A, n, a)
bench::mark(
  hydrorecipes:::exp_2(t, A, a),
  hydrorecipes:::exp_2_old(t, A, a),
  hydrorecipes:::gamma_3(t, A, n, a),
  hydrorecipes:::gamma_3_old(t, A, n, a),
  hydrorecipes:::gamma_3_old2(t, A, n, a),
  check = TRUE
)

bench::mark(
  hydrorecipes:::gamma_3(t, A, n, a),
  hydrorecipes:::gamma_3_old(t, A, n, a),
  hydrorecipes:::gamma_3_old2(t, A, n, a),
  check = TRUE
)


bench::mark(
  # hydrorecipes:::exp_2(t, A, a),
  # hydrorecipes:::exp_2_test(A, a),
  # hydrorecipes:::exp_2_old(t, A, a),
  # hydrorecipes:::gamma_3(t, A, n, a),
  hydrorecipes:::gamma_3(0:7200000, A, n, a),
  hydrorecipes:::gamma_3_test(A, n, a),
  hydrorecipes:::gamma_3_eigen(A, n, a),
  hydrorecipes:::exp_2(0:7200000, A, a),
  hydrorecipes:::exp_2_eigen(A, a),
  # hydrorecipes:::gamma_3_old(t, A, n, a),
  # hydrorecipes:::gamma_3_old2(t, A, n, a),
  check = TRUE
)

# x <- 1.0
#
# ig(rep(3.0, 1), 1.0)
#
#
# y <- rev(sort(abs(rnorm(1000000))))
# x <- 1/y
# bench::mark(
# hydrorecipes:::binary_search(x, y)
# )
#
# x <- rnorm(1e6)
# y <- 1.2
# bench::mark(
# hydrorecipes:::well_function_coefficient_rcpp(x,y),
# hydrorecipes:::well_function_coefficient_vec(x,y),
# check = FALSE
# )
#
# x <- abs(rnorm(1000000))
# bench::mark(
#             hydrorecipes:::ei_sp_vec(x),
#             hydrorecipes:::ei_vec(x),
#             hydrorecipes:::ei_eigen(x),
#             iterations = 1)
#
# n <- 1e5
# storativity = 1e-5
# radius = 50
# transmissivity = 1e-3
# leakage <- 100
# times <- as.numeric(1:n)
# flow_rates <- rep(1, n)#abs(rnorm(n))
# n_terms <- 12L
#
# bench::mark(
#
#   tmp <- hydrorecipes:::hantush_jacob(
#     times,
#     flow_rates,
#     radius,
#     storativity,
#     transmissivity,
#     leakage,
#     n_terms)[[1]],
#
#   tmp2 <- aquifer:::hantush_convolve(radius,
#                                      storativity,
#                                      transmissivity,
#                                      leakage,
#                                      times,
#                                      flow_rates,
#                                      n_terms),
#   check = FALSE
# )
#
# # x <- rnorm(100000)
# # bench::mark(hydrorecipes:::convolve_overlap_add(x,x),
# #             hydrorecipes:::convolve_filter(x,x,TRUE, TRUE), check = FALSE)
#
# plot(tmp2[1:1000], col = 'red', type = 'l', log = 'xy')
# points(tmp[1:1000], type = 'l', log = 'xy')
#
# x <- sort(abs(rnorm(500)))
# bench::mark(
#   # hydrorecipes:::hantush_well_vec(x, 0.01, 10),
#   # hydrorecipes:::hantush_well_e(x, 0.01, 10),
#   hydrorecipes:::hantush_well_rcpp(x, 0.01, 10),
#   aquifer:::hantush_well_parallel(x, 0.01, 10),
#   check = FALSE
# )
#
#
# bench::mark(
#   hydrorecipes:::hantush_jacob(x, rep(0.01, length(x)),
#                            radius = 10, 1e-6, 1e-3, 1, 20),
#   sapply(1:100, function(x) aquifer:::hantush_well_single(0.1, 0.01, 10)),
#   check = FALSE
# )
#
# x <- abs(rnorm(10000))
# bench::mark(
#   hydrorecipes:::bessel_eigen(x),
#   hydrorecipes:::bessel_rcpp(x),
#   besselK(x, 1),
#   Bessel::BesselK(x, 1)
# )
#
# library(expint)
# x <- abs(rnorm(10000))
# bench::mark(hydrorecipes:::gis(x, 0.0),
#             hydrorecipes:::gis2(x),
#             check = FALSE)
#
# n <- 100
# time <- seq(1, n, 1)
# flow_rate <- rep(0.001, n)
# thickness <- 10
# radius <- 5
# specific_storage_1 <- 1e-5
# hydraulic_conductivity_1 <- 1e-5
# specific_storage_2 <- 5e-6
# hydraulic_conductivity_2 <- 5e-6
# diffusivity_1 <- hydraulic_conductivity_1 / specific_storage_1
# diffusivity_2 <- hydraulic_conductivity_2 / specific_storage_2
# bench::mark(
#
# grf_1 <- hydrorecipes:::grf_time(radius,
#                              specific_storage_1,
#                              hydraulic_conductivity_1,
#                              thickness,
#                              time,
#                              flow_rate,
#                              flow_dimension = 2)
# )
# plot(grf_1[[1]], type = 'l', log = 'xy')
# points(grf_1[[1]], type = 'l', col = 'red')
# points(grf_1[[1]], type = 'l', col = 'blue')
# specific_storage <- 0.5e-5
# hydraulic_conductivity <- 0.5e-5
# grf_2 <- grf_time(radius,
#                   specific_storage_2,
#                   hydraulic_conductivity_2,
#                   thickness,
#                   time,
#                   flow_rate,
#                   flow_dimension = 2)
#
# dat <- data.table(time, grf_1, grf_2)
# fit_1 <- lm(grf_1~log(time), tail(dat, 4000))
# fit_2 <- lm(grf_2~log(time), tail(dat, 4000))
# summary(fit_1)
# summary(fit_2)
# plot(grf_2, log = 'x', type = 'l', xlab = "Elapsed time", ylab = "drawdown")
# abline(h = 0, col = 'grey')
# points(grf_1, type = 'l', col = 'red')
#
# dat[, pred_1 := predict(fit_1, dat)]
# dat[, pred_2 := predict(fit_2, dat)]
#
# points(pred_1~time, dat, type = 'l', col = 'red', lty = 2)
# points(pred_2~time, dat, type = 'l', lty = 2)
#
# hydrorecipes:::eig(1,1)
# library(expint)
# gammainc(1,1)
# gammainc(-1,1)
# hydrorecipes:::eig(1,1)
# hydrorecipes:::eig(1,0)


xy <- expand.grid(1:100, 1:100)
bench::mark(

  # Eigen::VectorXd distance_x,
  # Eigen::VectorXd distance_y,
  # Eigen::VectorXd output_times,
  # Eigen::VectorXd start_times,
  # Eigen::VectorXd flow_rates,
  # Eigen::VectorXd well_x,
  # Eigen::VectorXd well_y,
  # double storativity,
  # double transmissivity_x,
  # double transmissivity_y,
  # double thickness)

  a <- hydrorecipes:::grid_pumping_regimes(
    distance_x = xy[,1],
    distance_y = xy[,2],
    output_times = seq(1, 100, 1),
    start_times = sort(runif(min = 0, max = 90, 100)),
    flow_rates = rnorm(100),
    well_x = rep(c(200, 400, 500, 20, 10, 50, 7, 700, 800, 177), each = 10),
    well_y = rep(500, 100),
    storativity = 1e-6,
    transmissivity_x = 1e-4,
    transmissivity_y = 1e-5,
    thickness = 1.0)
)



*/
