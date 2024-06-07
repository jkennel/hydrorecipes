
#include "hydrorecipes.h"


// anisotropy value greater than 0
// major_axis_angle
// [[Rcpp::export]]
Eigen::MatrixXd coordinate_transform(Eigen::Map<Eigen::MatrixXd> coords,
                                     const double anisotropy,
                                     const double major_axis_angle) {

  if (anisotropy < 1.0) {
    Rcpp::stop("anisotropy must be greater than or equal to 1");
  }

  Eigen::Matrix2d rot;
  rot <<  cos(major_axis_angle), sin(major_axis_angle),
         -sin(major_axis_angle), cos(major_axis_angle);


  Eigen::Vector2d aniso;
  aniso << 1.0, 1.0 / anisotropy;

  return (coords * (rot * aniso.asDiagonal()));

}

// [[Rcpp::export]]
Eigen::MatrixXd coordinate_rotate(Eigen::Map<Eigen::MatrixXd> coords,
                                  const double major_axis_angle) {


  Eigen::Matrix2d rot;
  rot <<  cos(major_axis_angle), sin(major_axis_angle),
          -sin(major_axis_angle), cos(major_axis_angle);



  return (coords * rot);

}

/*** R
m <- matrix(c(10,10), ncol = 2)
bench::mark(
hydrorecipes:::coordinate_transform(m, 10.0, pi/2),
geoR::coords.aniso(m, c(pi/2, 10.0))
)
*/
