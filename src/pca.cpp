#include "frecipes.h"

// cite https://github.com/AEBilgrau/correlateR/

// // [[Rcpp::export]]
// Eigen::MatrixXd scale_eigen(Eigen::Map<Eigen::MatrixXd> & X) {
//
//   const double df = X.rows() - 1; // Subtract 1 by default
//
//   X.rowwise() -= X.colwise().mean();  // Centering
//
//   Eigen::MatrixXd cor = X.transpose() * X / df;   // The covariance matrix
//
//   // Get 1 over the standard deviations
//   Eigen::RowVectorXd inv_sds = cor.diagonal().array().sqrt().inverse();
//
//   return(X.array().rowwise() * inv_sds.array());
//
// }
//
//
// // [[Rcpp::export]]
// Rcpp::List scale_list_eigen(const Rcpp::List x,
//                             bool center = true,
//                             bool scale = true) {
//
//   size_t nc = x.size();
//   Eigen::VectorXd v_tmp = x[0];
//   size_t nr = v_tmp.size();
//   double n = double(nr) - 1.0;
//
//   Rcpp::List centered;
//
//   Eigen::VectorXd m_vector(nc);
//   Eigen::VectorXd sd_vector(nc);
//
//   for (int i = 0; i < nc; ++i) {
//
//     v_tmp = x[i];
//     m_vector(i) = v_tmp.mean();
//     v_tmp = v_tmp.array() - m_vector(i);
//     sd_vector(i) = sqrt(v_tmp.array().square().sum() / n);
//
//     centered.push_back(v_tmp.array() / sd_vector(i));
//
//   }
//
//   return Rcpp::List::create(
//     Rcpp::Named("x_scaled") = centered,
//     Rcpp::Named("center") = m_vector,
//     Rcpp::Named("scale") = sd_vector
//   );
//
// }
//
// // [[Rcpp::export]]
// Rcpp::List scale_param(const Rcpp::List x) {
//
//   size_t nc = x.size();
//   int nr;
//   double n;
//
//   Eigen::VectorXd m_vector(nc);
//   Eigen::VectorXd sd_vector(nc);
//
//   for (int i = 0; i < nc; ++i) {
//
//     Eigen::VectorXd v_tmp = x[i];
//     n = double(nr) - 1.0;
//     nr = v_tmp.size();
//     m_vector(i) = v_tmp.mean();
//     v_tmp = v_tmp.array() - m_vector(i);
//     sd_vector(i) = sqrt(v_tmp.array().square().sum() / n);
//
//   }
//
//   return Rcpp::List::create(
//     Rcpp::Named("center") = m_vector,
//     Rcpp::Named("scale") = sd_vector
//   );
//
// }
//

// [[Rcpp::export]]
Rcpp::List scale_list_param(const Rcpp::List x,
                            const NumericVector center,
                            const NumericVector scale) {

  unsigned int nc = center.size();
  Rcpp::NumericVector v_tmp;

  Rcpp::List centered(nc);

  for (unsigned int i = 0; i < nc; ++i) {
    v_tmp = Rcpp::clone(Rcpp::as<Rcpp::NumericVector>(x[i]));

    if (center(i) == 0.0){
      v_tmp = v_tmp * scale(i);
    } else if (scale(i) == 1.0){
      v_tmp = v_tmp - center(i);
    } else if (scale(i) == 1.0 & center(i) == 0.0) {
      v_tmp = (v_tmp - center(i)) * scale(i);
    }

    centered[i] = v_tmp;
  }

  return centered;

}

// [[Rcpp::export]]
Rcpp::List scale_list_param_std(const Rcpp::List x,
                            const NumericVector center,
                            const NumericVector scale) {

  unsigned int nc = center.size();
  std::vector<double> v_tmp;

  Rcpp::List centered(nc);

  for (unsigned int i = 0; i < nc; ++i) {
    v_tmp = Rcpp::as<std::vector<double> >(x[i]);

    if (scale(i) == 1.0 & center(i) == 0.0) {
    } else if (center(i) == 0.0){
      for (auto &out : v_tmp)
        out *= scale(i);
    } else if (scale(i) == 1.0){
      for (auto &out : v_tmp)
        out -= center(i);
    } else {
      for (auto &out : v_tmp)
        out = (out - center(i)) * scale(i);
    }

    centered[i] = v_tmp;
  }

  return centered;

}

// [[Rcpp::export]]
Rcpp::List scale_list_param_eigen(Rcpp::List x,
                                  Eigen::VectorXd center,
                                  Eigen::VectorXd scale) {

  size_t nc = center.size();
  Eigen::VectorXd v_tmp;

  Rcpp::List centered;

  for (size_t i = 0; i < nc; ++i) {
    v_tmp = x[i];

    if (scale(i) == 1.0 & center(i) == 0.0) {
    } else if (center(i) == 0.0){
      v_tmp = v_tmp.array() * scale(i);
    } else if (scale(i) == 1.0){
      v_tmp = v_tmp.array() - center(i);
    } else {
      v_tmp = (v_tmp.array() - center(i)) * scale(i);
    }

    centered.push_back(v_tmp);
  }

  return centered;

}

// [[Rcpp::export]]
Eigen::MatrixXd cor_list_eigen(Rcpp::List x,
                               Eigen::VectorXd center,
                               Eigen::VectorXd scale) {

  int nc = center.size();
  Rcpp::NumericVector tmp = x[0];
  int nr = tmp.size();
  double n = double(nr) - 1.0;

  Eigen::MatrixXd centered(nr, nc);
  Eigen::VectorXd v_tmp(nr);

  for (int i = 0; i < nc; ++i) {
    v_tmp = x[i];
    centered.col(i) = (v_tmp.array() - center(i)) / scale(i);
  }

  Eigen::MatrixXd cor = ((centered.adjoint() * centered) / n);

  return(cor);

}

// [[Rcpp::export]]
Rcpp::List pca_list_eigen(Rcpp::List x,
                          Eigen::RowVectorXd center,
                          Eigen::RowVectorXd scale,
                          bool prep = true) {

  int nc = center.size();
  Rcpp::NumericVector tmp = x[0];
  int nr = tmp.size();
  double n = double(nr) - 1.0;

  Eigen::MatrixXd centered(nr, nc);
  Eigen::VectorXd v_tmp(nr);

  for (int i = 0; i < nc; ++i) {
    v_tmp = x[i];
    centered.col(i) = (v_tmp.array() - center(i)) / scale(i);
  }

  Eigen::MatrixXd cor = ((centered.adjoint() * centered) / n);

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> pca(cor);

  Eigen::MatrixXd e_vectors = pca.eigenvectors().rowwise().reverse();

  if (prep) {
    return Rcpp::List::create(Rcpp::Named("rotation") = e_vectors);
  }

  Eigen::VectorXd e_values = pca.eigenvalues().reverse().array().sqrt();

  return Rcpp::List::create(
    Rcpp::Named("sdev") = e_values,
    Rcpp::Named("rotation") = e_vectors,
    Rcpp::Named("x") = centered * e_vectors,
    Rcpp::Named("cor") = cor
  );

}

// [[Rcpp::export]]
Eigen::MatrixXd scale_list_matrix_eigen(const Rcpp::List x,
                                        bool center = true,
                                        bool scale = true) {

  size_t nc = x.size();
  Eigen::VectorXd v_tmp = x[0];
  size_t nr = v_tmp.size();
  double n = double(nr) - 1.0;


  Eigen::VectorXd m_vector(nc);
  Eigen::VectorXd sd_vector(nc);

  Eigen::MatrixXd centered(nr, nc);
  for (int i = 0; i < nc; ++i) {

    v_tmp = x[i];
    m_vector(i) = v_tmp.mean();
    v_tmp = v_tmp.array() - m_vector(i);
    sd_vector(i) = sqrt(v_tmp.array().square().sum() / n);

    centered.col(i) = (v_tmp.array() / sd_vector(i));

  }

  return centered;


}

// [[Rcpp::export]]
Eigen::MatrixXd pca_list_rotation_eigen(Rcpp::List x,
                                        Eigen::VectorXd center,
                                        Eigen::VectorXd scale,
                                        int n_comp) {

  int nc = x.size();
  Rcpp::NumericVector tmp = x[0];
  int nr = tmp.size();
  double n = double(nr) - 1.0;

  Eigen::MatrixXd centered(nr, nc);
  Eigen::VectorXd v_tmp(nr);

  for (int i = 0; i < nc; ++i) {

    v_tmp = x[i];

    if (scale(i) == 1.0 & center(i) == 0.0) {
      centered.col(i) = v_tmp;
    } else if (center(i) == 0.0){
      centered.col(i) = v_tmp.array() / scale(i);
    } else if (scale(i) == 1.0){
      centered.col(i) = v_tmp.array() - center(i);
    } else {
      centered.col(i) = (v_tmp.array() - center(i)) / scale(i);
    }

  }

  Eigen::MatrixXd cor = (centered.adjoint() * centered) / n;

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> pca(cor);

  Eigen::MatrixXd e_vectors = pca.eigenvectors().rightCols(n_comp).rowwise().reverse();

  return e_vectors;

}


// [[Rcpp::export]]
Rcpp::List pca_eigen(const Rcpp::List x,
                     bool center = true,
                     bool scale = true) {

  size_t nc = x.size();
  Eigen::VectorXd v_tmp = x[0];
  size_t nr = v_tmp.size();
  double n = double(nr) - 1.0;


  Eigen::RowVectorXd m_vector(nc);
  Eigen::RowVectorXd sd_vector(nc);

  Eigen::MatrixXd centered(nr, nc);

  for (int i = 0; i < nc; ++i) {

    v_tmp = x[i];
    m_vector(i) = v_tmp.mean();
    v_tmp = v_tmp.array() - m_vector(i);
    sd_vector(i) = sqrt(v_tmp.array().square().sum() / n);

    centered.col(i) = (v_tmp.array() / sd_vector(i));

  }

  const Eigen::MatrixXd cor = ((centered.adjoint() * centered) / n);

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> pca(cor);

  Eigen::VectorXd e_values = pca.eigenvalues().reverse();
  const Eigen::MatrixXd e_vectors = pca.eigenvectors().rowwise().reverse();


  return Rcpp::List::create(
    Rcpp::Named("sdev") = e_values.array().sqrt(),
    Rcpp::Named("rotation") = e_vectors,
    Rcpp::Named("center") = m_vector,
    Rcpp::Named("scale") = sd_vector,
    Rcpp::Named("x") = centered * e_vectors
  );

}

// [[Rcpp::export]]
Eigen::MatrixXd cor_eigen(Eigen::Map<Eigen::MatrixXd> & X) {

  const double df = X.rows() - 1; // Subtract 1 by default

  X.rowwise() -= X.colwise().mean();  // Centering

  Eigen::MatrixXd cor = X.transpose() * X / df;   // The covariance matrix

  // Get 1 over the standard deviations
  Eigen::RowVectorXd inv_sds = cor.diagonal().array().sqrt().inverse();

  return(cor.array().rowwise() * inv_sds.array());

}


// // [[Rcpp::export]]
// Eigen::MatrixXd cor_list_eigen(Rcpp::List & x) {
//
//
//   Eigen::MatrixXd m_scale = scale_list_matrix_eigen(x, true, true);
//   double n = double(m_scale.rows()) - 1.0;
//
//   return(m_scale.transpose() * m_scale / n);
//
// }


// // [[Rcpp::export]]
// Eigen::MatrixXd scale_with_param_eigen(const Rcpp::List x,
//                                        Eigen::VectorXd center,
//                                        Eigen::VectorXd scale) {
//
//   int nc = x.size();
//   Eigen::VectorXd v_tmp = x[0];
//   int nr = v_tmp.size();
//   double n = double(nr) - 1.0;
//
//   Eigen::MatrixXd centered(nr, nc);
//
//   Eigen::VectorXd m_vector(nc);
//   Eigen::VectorXd sd_vector(nc);
//
//   for (int i = 0; i < nc; ++i) {
//     v_tmp = x[i];
//     centered.col(i) = (v_tmp.array() - center(i)) / scale(i);
//   }
//
//   return centered;
//
// }




// [[Rcpp::export]]
Rcpp::List pca(Eigen::Map<Eigen::MatrixXd> x,
               const bool center = true,
               const bool scale = true) {

  const double n = double(x.rows()) - 1.0;
  const int nc = x.cols();
  const int nr = x.rows();

  Eigen::VectorXd m_vector = Eigen::VectorXd::Zero(nc);
  Eigen::VectorXd sd_vector = Eigen::VectorXd::Zero(nc);
  Eigen::MatrixXd centered = x;

  for (int i = 0; i < nc; ++i) {

    if (center) {
      m_vector(i) = x.col(i).mean();
      centered.col(i) = x.col(i).array() - m_vector(i);
    }

    if (scale) {
      sd_vector(i) = centered.col(i).array().square().sum() / n;
      centered.col(i) = centered.col(i).array() / sqrt(sd_vector(i));
    }
  }

  Eigen::MatrixXd cor = ((centered.adjoint() * centered) / n);

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> pca(cor);

  Eigen::VectorXd e_values = pca.eigenvalues().reverse();
  Eigen::MatrixXd e_vectors = pca.eigenvectors().rowwise().reverse();


  return Rcpp::List::create(
    Rcpp::Named("sdev") = e_values.array().sqrt(),
    Rcpp::Named("rotation") = e_vectors,
    Rcpp::Named("center") = m_vector,
    Rcpp::Named("scale") = sd_vector.array().sqrt(),
    Rcpp::Named("x") = centered * e_vectors
  );

}




// [[Rcpp::export]]
Rcpp::List pca_with_params(Eigen::Map<Eigen::MatrixXd> x,
                           const Eigen::RowVectorXd center,
                           const Eigen::RowVectorXd scale) {

  const double n = double(x.rows()) - 1.0;
  const int nc = x.cols();
  const int nr = x.rows();

  Eigen::MatrixXd centered = x;

  for (int i = 0; i < nc; ++i) {
    centered.col(i) = (x.col(i).array() - center(i)).array() / scale(i);
  }

  Eigen::MatrixXd cor = ((centered.adjoint() * centered) / n);

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> pca(cor);

  Eigen::VectorXd e_values = pca.eigenvalues().reverse();
  Eigen::MatrixXd e_vectors = pca.eigenvectors().rowwise().reverse();


  return Rcpp::List::create(Rcpp::Named("rotation") = e_vectors,
                            Rcpp::Named("sdev") = e_values.array().sqrt(),
                            Rcpp::Named("center") = center,
                            Rcpp::Named("scale") = scale,
                            Rcpp::Named("x") = centered * e_vectors

  );
}

// [[Rcpp::export]]
Rcpp::List pca_list_with_params(Rcpp::List x,
                                Eigen::RowVectorXd center,
                                Eigen::RowVectorXd scale) {

  int nc = x.size();
  Rcpp::NumericVector tmp = x[0];
  int nr = tmp.size();
  double n = double(nr) - 1.0;

  Eigen::MatrixXd centered(nr, nc);
  Eigen::VectorXd v_tmp(nr);

  for (int i = 0; i < nc; ++i) {
    v_tmp = x[i];

    centered.col(i) = (v_tmp.array() - center(i)).array() / scale(i);
  }

  Eigen::MatrixXd cor = ((centered.adjoint() * centered) / n);

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> pca(cor);

  Eigen::VectorXd e_values = pca.eigenvalues().reverse();
  Eigen::MatrixXd e_vectors = pca.eigenvectors().rowwise().reverse();


  return Rcpp::List::create(Rcpp::Named("rotation") = e_vectors,
                            Rcpp::Named("sdev") = e_values.array().sqrt(),
                            Rcpp::Named("center") = center,
                            Rcpp::Named("scale") = scale,
                            Rcpp::Named("x") = centered * e_vectors

  );
}



// // [[Rcpp::export]]
// Eigen::MatrixXd corEigen(Eigen::Map<Eigen::MatrixXd> & X) {
//
//   // Handle degenerate cases
//   if (X.rows() == 0 && X.cols() > 0) {
//     return Eigen::MatrixXd::Constant(X.cols(), X.cols(),
//                                      Rcpp::NumericVector::get_na());
//   }
//
//   // Computing degrees of freedom
//   // n - 1 is the unbiased estimate whereas n is the MLE
//   const int df = X.rows() - 1; // Subtract 1 by default
//
//   X.rowwise() -= X.colwise().mean();  // Centering
//
//   Eigen::MatrixXd cor = X.transpose() * X / df;   // The covariance matrix
//
//   // Get 1 over the standard deviations
//   Eigen::VectorXd inv_sds = cor.diagonal().array().sqrt().inverse();
//
//   // Scale the covariance matrix
//   cor = cor.cwiseProduct(inv_sds * inv_sds.transpose());
//
//   return cor;
// }
//
// // [[Rcpp::export]]
// Eigen::MatrixXd corEigen2(Eigen::Map<Eigen::MatrixXd> & X) {
//
//   // Handle degenerate cases
//   if (X.rows() == 0 && X.cols() > 0) {
//     return Eigen::MatrixXd::Constant(X.cols(), X.cols(),
//                                      Rcpp::NumericVector::get_na());
//   }
//
//   // Computing degrees of freedom
//   // n - 1 is the unbiased estimate whereas n is the MLE
//   const int df = X.rows() - 1; // Subtract 1 by default
//
//   X.rowwise() -= X.colwise().mean();  // Centering
//
//   Eigen::MatrixXd cor = X.adjoint() * X / df;   // The covariance matrix
//
//   // Get 1 over the standard deviations
//   Eigen::VectorXd inv_sds = cor.diagonal().array().sqrt().inverse();
//
//   // Scale the covariance matrix
//   cor = cor.cwiseProduct(inv_sds * inv_sds.transpose());
//
//   return cor;
// }


// // [[Rcpp::export]]
// Eigen::MatrixXd svd_eigen_jac(Eigen::MatrixXd C)
// {
//
//     Eigen::JacobiSVD<Eigen::MatrixXd> svd(C, ComputeThinU | ComputeThinV);
//     MatrixXd Cp = svd.matrixV();
//
//    return(Cp);
// }
//
// // [[Rcpp::export]]
// Eigen::MatrixXd svd_eigen_bdc(Eigen::MatrixXd C)
// {
//
//     Eigen::BDCSVD<Eigen::MatrixXd> svd(C, Eigen::ComputeThinU | Eigen::ComputeThinV);
//     MatrixXd Cp = svd.matrixV().leftCols(5);
//
//    return(Cp);
// }

/*** R
library(frecipes)
library(collapse)
nc <- 100
m <- matrix(rnorm(1e7), ncol = nc)
l <- unclass(qDF(m))
center <- fmean(l)
scale <- fsd(l)
bench::mark(frecipes:::scale_list_param_eigen(l, center, scale),
            frecipes:::scale_list_param(l, center, scale),
            frecipes:::scale_list_param_std(l, center, scale)
            )

nc <- 100
m <- matrix(rnorm(2e6), ncol = nc)
l <- unclass(qDF(m))
center <- fmean(m)
scale <- fsd(m)

bench::mark(
a <- svd(m)$u,
aa <- corpcor::fast.svd(m),
b <- frecipes:::svd_eigen_jac(m),
d <- frecipes:::svd_eigen_bdc(m),
aaaa <- frecipes:::pca_list_eigen(l, rep(0, nc), rep(1, nc), FALSE),
check = FALSE
)


library(frecipes)
library(collapse)
m <- as.matrix(USArrests)
m <- matrix(rnorm(2e6), ncol = 10)
l <- unclass(qDF(m))
center <- fmean(m)
scale <- fsd(m)
bench::mark(
  a <- frecipes:::pca(m),
  aa <- frecipes:::pca_list(l),
  aaaa <- frecipes:::pca_list_rotation(l, center, scale),
  d <- frecipes:::pca_with_params(m, center, scale),
  dd <- frecipes:::pca_list_with_params(l, center, scale),
  # ddd <- frecipes:::pca_list_with_params(l, fmean(l), fsd(l)),
  b <- prcomp(m,
              retx = FALSE,
              center = TRUE,
              scale = TRUE,
              tol = NULL,
              rank = 2),
  # e <- prcomp(m,
  #             retx = TRUE,
  #             center = TRUE,
  #             scale = TRUE,
  #             tol = NULL),
  # f <- prcomp(m,
  #             retx = TRUE,
  #             center = TRUE,
  #             scale = TRUE,
  #             tol = NULL),
  # g <- princomp(m,
  #               cor = TRUE,
  #               fix_sign = FALSE),
  check = FALSE)
head(a[[1]])
head(as.matrix(b$rotation))
(a[[2]])
(b$sdev)

fmean(m)
a[[3]]

(fsd(m))
(a[[4]])

rec <- recipe(~., data = as.data.frame(m))
pca_trans <- rec %>%
  step_normalize(all_numeric()) %>%
  step_pca(all_numeric(), num_comp = 3)
pca_estimates <- prep(pca_trans, training = USArrests)
pca_data <- bake(pca_estimates, USArrests)

mdf <- as.data.frame(m)
bench::mark(
  # {rec <- recipe(~., data = mdf)  |>
  #   step_pca(all_numeric(), num_comp = 3) |>
  #   prep() |>
  #   bake(new_data = NULL)},
  a <- frecipes:::pca(m, center = FALSE, scale = FALSE),
  b <- prcomp(m,
              retx = TRUE,
              center = FALSE,
              scale = FALSE,
              tol = NULL),
  check = FALSE

)
m <- collapse::qM(dat)


m <- matrix(rnorm(2e6), ncol = 100)
l <- unclass(qDF(m))
c <- fmean(l)
s <- fsd(l)
bench::mark(

  a <- fscale(m),
  qM(fscale(l)),
  # fsd(l,stable.algo=FALSE),
  # fmean(l,stable.algo=FALSE),
  scale(m, TRUE, TRUE),
  d <- frecipes:::scale_list_to_matrix_2(l, fmean(l), fsd(l)),
  b <- frecipes:::scale_list_to_matrix(l),
  check = FALSE

)


{rec <- recipe(~., data = mdf)  |>
    step_pca(all_numeric(), num_comp = 3) |>
    prep()}


library(correlateR)
n <- 50
m <- matrix(rnorm(2e6), ncol = n)
l <- unclass(qDF(m))
bench::mark(
  # a <- fscale(m),
  # corEigen(m),
  # e <- frecipes:::cor_eigen(m),
  # ee <- frecipes:::pca_eigen(l),
  b <- fscale(l),
  d <- frecipes:::scale_list_eigen(l),
  d2 <- frecipes:::scale_list_param_eigen(l, d$center, d$scale),
  dd <- frecipes:::scale_list_matrix_eigen(l),
  ddd <- qM(frecipes:::scale_list_eigen(l)[[1]]),
  e <- frecipes:::cor_list_eigen(l),
  f <- frecipes:::pca_list_eigen(l, scale = TRUE, center = TRUE, prep = FALSE),
  ff <- frecipes:::pca_list_eigen(l, scale = TRUE, center = TRUE, prep = TRUE),
  fff <- qM(l),
  # eeee <- frecipes:::pca(m),
  # b <- frecipes:::scale_eigen(m),
  # f <- frecipes:::scale_list_matrix_eigen(l),
  g <- prcomp(m, center = TRUE, scale = TRUE),
  h <- cor(m),
  check = FALSE,
  relative = FALSE
)

bench::mark(
  a <- frecipes:::scale_param(l),
  b <- frecipes:::scale_list_eigen(l),
  d <- collapse::fsd(l),
  d <- collapse::fmean(l),
  e <- lapply(l,sd),
  f <- lapply(l,mean),
  check = FALSE
)

*/
