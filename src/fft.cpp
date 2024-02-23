#include "frecipes.h"

//==============================================================================
//' @title
 //' fft_matrix
 //'
 //' @description
 //' Do an FFT for each matrix column
 //'
 //' @param x the matrix that holds the series (numeric matrix)
 //' @param detrend remove the linear trend of the columns (boolean)
 //' @param demean remove the mean for each column (boolean)
 //' @param n_new the padded size (integer)
 //'
 //' @return A matrix with FFT results.
 //'
 //' @noRd
 //'
 // [[Rcpp::export]]
 Eigen::MatrixXcd fft_matrix(Eigen::MatrixXd x,
                             size_t n_new) {

   size_t n_row = x.rows();
   size_t n_col = x.cols();

   VectorXcd x_fft(n_new);
   VectorXd  x_padded(n_new);
   MatrixXcd out(n_new, n_col);

   Eigen::FFT<double> fft;

   for (size_t i = 0; i < n_col; ++i) {
     x_padded = pad_vector(x.col(i), n_row, n_new);
     fft.fwd(x_fft, x_padded, 0);
     out.col(i) = x_fft;
   }


   return(out);
 }
//==============================================================================




//******************************************************************************
// Convolution
//******************************************************************************
//==============================================================================
//' @title
 //' convolve_vec
 //'
 //' @description
 //' Circular convolution of two vectors having the same length
 //'
 //' @param x the vector that holds the series (numeric vector)
 //' @param y the vector to convolve with x (numeric vector)
 //'
 //'
 //' @return numeric vector that is the circular convolution of two vectors
 //'
 //'
 //' @noRd
 //'
 // [[Rcpp::export]]
 Eigen::VectorXd convolve_vec(Eigen::VectorXd x,
                              Eigen::VectorXd y) {

   Eigen::FFT<double> fft;
   size_t n_x = x.size();
   size_t n_y = y.size();

   if (n_y != n_x) {
     Rcpp::stop("convolve_vec: the lengths of x and y should be the same");
   }

   VectorXcd fft_x(n_x);
   VectorXcd fft_y(n_x);
   VectorXd z(n_x);

   fft.fwd(fft_x, x);
   fft.fwd(fft_y, y);

   fft_x = fft_x.array() * fft_y.array();
   fft.inv(z, fft_x);

   return(z);
 }
//==============================================================================

//==============================================================================
//' @title
 //' convolve_filter
 //'
 //' @description
 //' convolution of vector with matrix
 //'
 //' @param x vector to convolve with y (numeric vector)
 //' @param y numeric matrix to convolve with x (column by column convolution)
 //'  (numeric matrix)
 //' @param remove_partial keep the end values or fill with NA (boolean)
 //' @param reverse should x be reversed before convolution (boolean)
 //'
 //' @return numeric matrix of convolved values
 //'
 //' @export
 //'
 //' @importFrom Rcpp sourceCpp
 //' @importFrom stats nextn
 //' @importFrom stats convolve
 //' @importFrom stats spec.pgram
 //'
 //' @examples
 //' a <- convolve_filter(x = 1:100,
 //'                      y = c(1:10, rep(0, 90)),
 //'                      remove_partial = FALSE,
 //'                      reverse = TRUE)
 //'
 //' b <- stats::convolve(1:100, rev(1:10), type = 'filter')
 //'
 // [[Rcpp::export]]
 Eigen::VectorXd convolve_filter(const Eigen::VectorXd& x,
                                 const Eigen::VectorXd& y,
                                 const bool remove_partial,
                                 const bool reverse) {

   size_t n_x = x.size();
   size_t n_y = y.size();

   Eigen::FFT<double> fft;
   size_t n_new = next_n_eigen(n_x + n_y - 1);

   // Temporary vectors
   VectorXd x_dbl = VectorXd::Zero(n_new);
   VectorXd out_dbl = VectorXd::Zero(n_new);

   VectorXcd fft_x(n_new);
   VectorXcd fft_y(n_new);

   // Output
   VectorXd out(n_y);

   if (reverse) {
     x_dbl.tail(n_x) = x.reverse();
   } else {
     x_dbl.tail(n_x) = x;
   }

   // do fft for x
   fft.fwd(fft_x, x_dbl);

   // do fft for y and convolution
   out_dbl.head(n_y) = y;

   fft.fwd(fft_y, out_dbl);

   fft_y = fft_x.array() * fft_y.conjugate().array();

   fft.inv(out_dbl, fft_y);

   if (reverse) {
     out = out_dbl.tail(n_x).reverse();
   } else {
     out = out_dbl.head(n_x);
   }

   if (remove_partial) {
     out.head(n_y - 1).setConstant(NA_REAL);
   }

   return(out);

 }
//==============================================================================


//==============================================================================
//' @title
 //' convolve_overlap_add
 //'
 //' @description
 //' Multiply a transfer function with a real input and take the inverse FFT.
 //'
 //' @param x the vector that holds the series (numeric vector)
 //' @param y the kernel to convolve with x (complex numeric vector)
 //'
 //' @return the linear convolution of two vectors
 //'
 //' @noRd
 //'
 // [[Rcpp::export]]
 Eigen::VectorXd convolve_overlap_add(Eigen::VectorXd& x,
                                      Eigen::VectorXd& y) {


   Eigen::FFT<double> fft;
   size_t n_x = x.size();
   size_t n_y = y.size();
   size_t n_pad = next_n_eigen(n_y * 4);
   if (n_y > n_x) Rcpp::stop("n_y cannot be larger than n_x");
   if (n_pad > n_x) n_pad = n_x;

   // size_t n_times = n_x / n_y - 1;
   Eigen::VectorXd x_sub = Eigen::VectorXd::Zero(n_pad);
   Eigen::VectorXd y_sub = pad_vector(y.reverse(), n_y, n_pad);

   VectorXcd fft_y(n_pad);
   fft.fwd(fft_y, y_sub);
   VectorXcd fft_x(n_pad);
   VectorXd z(n_pad);
   VectorXd out = Eigen::VectorXd::Zero(n_x + n_pad);


   size_t i = 0;
   size_t x_len = n_pad - n_y + 1;

   while (i < n_x) {
     x_len = std::min(x_len, (n_x - i));

     x_sub.head(x_len) = x.segment(i, x_len);

     fft.fwd(fft_x, x_sub);
     fft_x = fft_x.array() * fft_y.array();
     fft.inv(z, fft_x);
     out.segment(i, n_pad) += z;

     i += x_len;

   }

   out.head(n_y - 1).setConstant(NA_REAL);

   // out.conservativeResize(n_x);

   return(out.head(n_x));
 }
//==============================================================================



//==============================================================================
//' @title
 //' convolve_overlap_save
 //'
 //' @description
 //' Multiply a transfer function with a real input and take the inverse FFT.
 //'
 //' @param x the vector that holds the series (numeric vector)
 //' @param y the kernel to convolve with x (complex numeric vector)
 //' @param align right (0), center (1), or left (2) alignment
 //'
 //'
 //' @return the linear convolution of two vectors
 //'
 //' @noRd
 //'
 // [[Rcpp::export]]
 Eigen::VectorXd convolve_overlap_save(Eigen::VectorXd x,
                                       Eigen::VectorXd y,
                                       int align) {

   Eigen::FFT<double> fft;
   size_t n_x = x.size();
   size_t n_y = y.size();

   // need a better way of choosing this
   size_t n_pad = next_n_eigen(n_y * 4);
   if (n_y > n_x) Rcpp::stop("n_y cannot be larger than n_x");
   if (n_pad > n_x) n_pad = n_x;


   Eigen::VectorXd x_sub = Eigen::VectorXd::Zero(n_pad);
   Eigen::VectorXd y_sub = pad_vector(y.reverse(), n_y, n_pad);

   VectorXcd fft_y(n_pad);
   fft.fwd(fft_y, y_sub);
   VectorXcd fft_x(n_pad);
   VectorXd z(n_pad);
   VectorXd out = VectorXd::Zero(n_x);


   size_t i = n_x - n_pad;
   size_t x_len = n_pad;
   size_t fin_size = x_len - n_y;
   size_t n_align = n_y; // right alignment

   if (align == 1) {
     n_align = (int)(n_y / 2);
   }
   if (align == 2) {
     n_align = 1;
   }

   while (i >= 0) {

     fft.fwd(fft_x, x.segment(i, x_len));
     fft_x = fft_x.array() * fft_y.array();
     fft.inv(z, fft_x);

     out.segment(i + n_align, fin_size + 1) = z.tail(fin_size + 1);


     if (i == 0) {
       break;
     } else if (fin_size > i) {
       i = 0;
     } else {
       i -= fin_size;
     }


   }

   if (align == 0) {
     out.head(n_y - 1).setConstant(NA_REAL);
   }
   if (align == 1) {
     out.head(n_align).setConstant(NA_REAL);
     out.tail(n_align).setConstant(NA_REAL);
   }
   if (align == 2) {
     out.tail(n_y - 1).setConstant(NA_REAL);
   }

   // out.conservativeResize(n_x);

   return(out);
 }
//==============================================================================

// [[Rcpp::export]]
Eigen::VectorXd shift_eigen(Eigen::VectorXd x, int n)
{

  std::vector<double> x_vec(x.data(), x.data() + x.size());


  std::rotate(x_vec.begin(),
              x_vec.end() - n, // this will be the new first element
              x_vec.end());

  x = Eigen::Map<Eigen::VectorXd>(x_vec.data(), x.size());

  return(x);
}

//==============================================================================
//' @title
 //' convolve_overlap_save_list
 //'
 //' @description
 //' Multiply a transfer function with a real input and take the inverse FFT.
 //'
 //' @param x the vector that holds the series (numeric vector)
 //' @param y the list of kernels to convolve with x
 //' @param align right (0), center (1), or left (2) alignment
 //'
 //' @return the linear convolution of two vectors
 //'
 //' @noRd
 //'
 // [[Rcpp::export]]
 List convolve_overlap_save_list(Eigen::VectorXd& x,
                                 List y,
                                 int align) {

   size_t n_y = y.size();
   size_t n_x = x.size();
   Eigen::VectorXd out(n_x);
   Rcpp::List out_list(n_y);

   for (size_t j = 0; j < n_y; ++j) {
     out = convolve_overlap_save(x, y[j], align);
     out_list[j] = out;
   }

   return(out_list);
 }


//==============================================================================
//' @title
 //' convolve_tf
 //'
 //' @description
 //' Multiply a transfer function with a real input and take the inverse FFT.
 //'
 //' @param x the vector that holds the series (numeric vector)
 //' @param y the transfer function to multiply with x (complex numeric vector)
 //'
 //'
 //' @return the circular convolution of two vectors
 //'
 //' @noRd
 //'
 // [[Rcpp::export]]
 Eigen::VectorXd convolve_tf(Eigen::VectorXd x,
                             Eigen::VectorXcd y) {

   Eigen::FFT<double> fft;
   size_t n_x = x.size();
   size_t n_y = y.size();

   if (n_y != n_x) {
     Rcpp::stop("convolve_tf: the lengths of x and y should be the same");
   }

   VectorXcd fft_x(n_x);
   VectorXd z(n_x);

   fft.fwd(fft_x, x);

   fft_x = fft_x.array() * y.array();
   fft.inv(z, fft_x);

   return(z);
 }
//==============================================================================


//==============================================================================
//' @title
 //' convolve_matrix
 //'
 //' @description
 //' convolution of vector with matrix
 //'
 //' @param x vector to convolve with y (numeric vector)
 //' @param y numeric matrix to convolve with x (column by column convolution)
 //'  (numeric matrix)
 //' @param remove_partial keep the end values or fill with NA (boolean)
 //' @param reverse should x be reversed before convolution (boolean)
 //'
 //' @return numeric matrix of convolved values
 //'
 //' @export
 //'
 //' @importFrom Rcpp sourceCpp
 //' @importFrom stats nextn
 //' @importFrom stats convolve
 //' @importFrom stats spec.pgram
 //'
 //' @examples
 //' a <- convolve_matrix(x = 1:100,
 //'                      y = as.matrix(1:10),
 //'                      remove_partial = FALSE,
 //'                      reverse = TRUE)
 //'
 //' b <- stats::convolve(1:100, rev(1:10), type = 'filter')
 //'
 // [[Rcpp::export]]
 Eigen::MatrixXd convolve_matrix(const Eigen::VectorXd& x,
                                 const Eigen::MatrixXd& y,
                                 const bool remove_partial,
                                 const bool reverse) {

   size_t n_row   = x.size();
   size_t n_row_y = y.rows();
   size_t n_col   = y.cols();

   if (n_row_y > n_row) {
     Rcpp::stop("convolve_matrix: y cannot have more rows than x");
   }

   Eigen::FFT<double> fft;
   size_t n_new = next_n_eigen(n_row + n_row_y - 1);

   // Temporary vectors
   VectorXd x_dbl(n_new);
   VectorXd out_dbl(n_new);

   VectorXcd fft_x(n_new);
   VectorXcd fft_y(n_new);

   // Output
   MatrixXd out(n_row, n_col);

   x_dbl.setZero();

   if (reverse) {
     x_dbl.tail(n_row) = x.reverse();
   } else {
     x_dbl.tail(n_row) = x;
   }

   // do fft for x
   fft.fwd(fft_x, x_dbl);


   // do fft for y and convolution
   for (size_t i = 0; i < n_col; ++i) {
     out_dbl.setZero();

     out_dbl.head(n_row_y) = y.col(i);

     fft.fwd(fft_y, out_dbl);

     fft_y = fft_x.array() * fft_y.conjugate().array();

     fft.inv(out_dbl, fft_y);

     if (reverse) {
       out.col(i) = out_dbl.tail(n_row).reverse();
     } else {
       out.col(i) = out_dbl.head(n_row);
     }
   }

   if (remove_partial) {
     out.topRows(n_row_y - 1).setConstant(NA_REAL);
   }

   return(out);

 }
//==============================================================================

//==============================================================================
//' @title
 //' convolve_matrix
 //'
 //' @description
 //' convolution of vector with matrix
 //'
 //' @param x vector to convolve with y (numeric vector)
 //' @param y numeric matrix to convolve with x (column by column convolution)
 //'  (numeric matrix)
 //' @param remove_partial keep the end values or fill with NA (boolean)
 //' @param reverse should x be reversed before convolution (boolean)
 //'
 //' @return numeric matrix of convolved values
 //'
 //' @export
 //'
 //' @importFrom Rcpp sourceCpp
 //' @importFrom stats nextn
 //' @importFrom stats convolve
 //' @importFrom stats spec.pgram
 //'
 //' @examples
 //' a <- convolve_matrix(x = 1:100,
 //'                      y = as.matrix(1:10),
 //'                      remove_partial = FALSE,
 //'                      reverse = TRUE)
 //'
 //' b <- stats::convolve(1:100, rev(1:10), type = 'filter')
 //'
 // [[Rcpp::export]]
 Rcpp::List convolve_list(const Eigen::VectorXd& x,
                          const List y,
                          const bool remove_partial,
                          const bool reverse) {

   size_t n = y.length();
   size_t n_x = x.size();

   Eigen::VectorXd yy = Rcpp::as<Eigen::Map<Eigen::VectorXd>>(y[0]);
   size_t n_y = yy.size();

   Eigen::FFT<double> fft;
   size_t n_new = next_n_eigen(n_x + n_y - 1);

   // Temporary vectors
   VectorXd x_dbl = VectorXd::Zero(n_new);
   VectorXd out_dbl = VectorXd::Zero(n_new);

   VectorXcd fft_x(n_new);
   VectorXcd fft_y(n_new);

   // Output
   Rcpp::List out;

   if (reverse) {
     x_dbl.tail(n_x) = x.reverse();
   } else {
     x_dbl.tail(n_x) = x;
   }

   // do fft for x
   fft.fwd(fft_x, x_dbl);

   // do fft for y and convolution

   for (size_t i = 0; i < n; ++i) {
     yy = Rcpp::as<Eigen::Map<Eigen::VectorXd>>(y[i]);

     out_dbl.head(n_y) = yy;

     fft.fwd(fft_y, out_dbl);

     fft_y = fft_x.array() * fft_y.conjugate().array();

     fft.inv(out_dbl, fft_y);

     if (reverse) {

       if (remove_partial) {
         out_dbl.tail(n_y - 1).setConstant(NA_REAL);
       }

       out.push_back(out_dbl.tail(n_x).reverse());

     } else {

       if (remove_partial) {
         out_dbl.head(n_y - 1).setConstant(NA_REAL);
       }

       out.push_back(out_dbl.head(n_x));

     }

   }


   return(out);

 }

//==============================================================================
// [[Rcpp::export]]
std::list<Eigen::VectorXd> convolve_list2(const Eigen::VectorXd& x,
                                          const std::list<Eigen::VectorXd> y,
                                          const bool remove_partial,
                                          const bool reverse) {

  size_t n = y.size();
  size_t n_x = x.size();

  Eigen::VectorXd yy = y.front();
  size_t n_y = yy.size();

  Eigen::FFT<double> fft;
  size_t n_new = next_n_eigen(n_x + n_y - 1);

  // Temporary vectors
  VectorXd x_dbl = VectorXd::Zero(n_new);
  VectorXd out_dbl = VectorXd::Zero(n_new);

  VectorXcd fft_x(n_new);
  VectorXcd fft_y(n_new);

  // Output
  std::list<Eigen::VectorXd> out;

  if (reverse) {
    x_dbl.tail(n_x) = x.reverse();
  } else {
    x_dbl.tail(n_x) = x;
  }

  // do fft for x
  fft.fwd(fft_x, x_dbl);

  // do fft for y and convolution
  for (auto yy : y) {
    // for (size_t i = y.cbegin(); i < y.cend(); ++i) {
    // yy = y[i];

    out_dbl.head(n_y) = yy;

    fft.fwd(fft_y, out_dbl);

    fft_y = fft_x.array() * fft_y.conjugate().array();

    fft.inv(out_dbl, fft_y);

    if (reverse) {

      if (remove_partial) {
        out_dbl.tail(n_y - 1).setConstant(NA_REAL);
      }

      out.push_back(out_dbl.tail(n_x).reverse());

    } else {

      if (remove_partial) {
        out_dbl.head(n_y - 1).setConstant(NA_REAL);
      }

      out.push_back(out_dbl.head(n_x));

    }

  }


  return(out);

}


//******************************************************************************
// Pgram
//******************************************************************************


//==============================================================================
//' @title
 //' multiply_ffts
 //'
 //' @description
 //' Multiply each column of a complex matrix with all the columns.
 //'
 //' @param x complex matrix to convolve with itself (complex numeric matrix)
 //' @param n_col number of columns in the original input series (integer)
 //' @param truncated skip the first row to decrease memory use? (boolean)
 //'
 //'
 //' @return pgram of input FFT values.
 //'
 //' @noRd
 //'
 // [[Rcpp::export]]
 Eigen::MatrixXcd multiply_ffts(Eigen::MatrixXcd& x) {

   size_t ind = 0;
   size_t n_col = x.cols();

   size_t n_col_pgram = double(n_col / 2.0) * (n_col + 1);

   size_t n_row = x.rows();
   MatrixXcd pgram_mat(n_row, n_col_pgram);

   // complete set
   for (size_t i = 0; i < n_col; ++i) {
     for (size_t j = i; j < n_col; ++j) {

       pgram_mat.col(ind) = x.col(i).array() *
         x.col(j).conjugate().array();

       ind += 1;

     }
   }

   return(pgram_mat);

 }
//==============================================================================


//==============================================================================
//' @title
 //' fill_lower_left
 //'
 //' @description
 //' Fill in the complex conjugate columns.
 //'
 //' @param x complex matrix of pgram values (complex matrix)
 //' @param n_col number of columns in the original input series (integer)
 //' @param start the first row index to begin on (boolean)
 //'
 //'
 //' @return Matrix with filled in complex conjugate columns.
 //'
 //' @noRd
 //'
 // [[Rcpp::export]]
 Eigen::MatrixXcd fill_lower_left(Eigen::MatrixXcd& x,
                                  size_t start) {
   size_t ind;
   size_t ind_inv;

   size_t n_col = sqrt(x.cols() * 2 - x.cols());

   // assign lower left with the complex conjugate
   for (size_t i = start; i < n_col - 1 + start; ++i){
     for (size_t j = i + 1; j < n_col; ++j){
       ind     = index_from_i_j(i - start, j, n_col);
       ind_inv = index_from_j_i(i, j - start, n_col);
       x.col(ind_inv) = x.col(ind).conjugate();
     }
   }

   return(x);
 }
//==============================================================================


//==============================================================================
//' @title
 //' spec_pgram
 //'
 //' @description
 //' Calculate the periodogram.  This method only keeps the columns necessary for
 //' the transfer function calculation. This method is based on `spec.pgram`.
 //'
 //' @inheritParams spec.pgram
 //'
 //'
 //' @return periodogram from an input matrix using a Fast Fourier Transform.
 //' Similar to `spec.pgram` but should be faster.
 //'
 //' @noRd
 //'
 // [[Rcpp::export]]
 Eigen::MatrixXcd spec_pgram(Eigen::MatrixXd& x,
                             const Eigen::VectorXi& spans,
                             bool detrend,
                             bool demean,
                             double taper) {

   if (spans.size() < 1) {
     Rcpp::stop("spec_pgram: spans must be length 1 or larger.");
   }

   // detrend or demean
   x = detrend_and_demean_matrix(x, detrend, demean);

   size_t n_row = x.rows();
   size_t n_new = next_n_eigen(n_row);

   std::complex<double> scale = 1.0 / n_row; // or n_new

   // taper vector
   ArrayXd taper_array = spec_taper(n_row, taper).array();


   // Do FFTs
   MatrixXcd x_fft_mat = fft_matrix(x.array().colwise() * taper_array,
                                    n_new);

   // do a check on small values
   // x_fft_mat = check_ffts(x_fft_mat, 1000);


   // calculate upper triangle
   MatrixXcd pgram_mat = multiply_ffts(x_fft_mat).array() * scale;


   // interpolate first value (do we want to do this?)
   pgram_mat.row(0) = 0.5 * (pgram_mat.row(1).array() + pgram_mat.row(n_new - 1).array());


   // kernel multiplication
   if (spans(0) > 1) {
     VectorXd kernel = modified_daniell(spans / 2);
     pgram_mat = kernel_apply(pgram_mat, kernel);
   }


   return(pgram_mat);
 }
//==============================================================================


//==============================================================================
//' @title
 //' spec_pgram
 //'
 //' @description
 //' Calculate the periodogram.  This method only keeps the columns necessary for
 //' the transfer function calculation. This method is based on `spec.pgram`.
 //'
 //' @inheritParams spec.pgram
 //'
 //'
 //' @return periodogram from an input matrix using a Fast Fourier Transform.
 //' Similar to `spec.pgram` but should be faster.
 //'
 //' @noRd
 //'
 // [[Rcpp::export]]
 Rcpp::List spec_pgram_list(Rcpp::List& x,
                            const Eigen::VectorXi& spans,
                            bool detrend,
                            bool demean,
                            double taper) {

   if (spans.size() < 1) {
     Rcpp::stop("spec_pgram: spans must be length 1 or larger.");
   }

   Eigen::VectorXd tmp = x[0];
   size_t ind = 0;
   size_t n_row = tmp.size();
   size_t n_new = next_n_eigen(n_row);
   size_t n = x.size();
   size_t n_col_pgram = double(n / 2.0) * (n + 1);
   Rcpp::List pgram(n_col_pgram);

   std::complex<double> scale = 1.0 / n_row; // or n_new

   // taper vector
   Eigen::ArrayXd taper_array = spec_taper(n_row, taper).array();

   // Do FFTs
   Rcpp::List x_fft_list(n);// = fft_list(x.array().colwise() * taper_array,n_new);

   Eigen::VectorXcd x_fft(n_new);
   Eigen::VectorXd  x_padded(n_new);
   Eigen::VectorXd  x_tmp(n_row);
   // MatrixXcd out(n_new, n_col);

   Eigen::FFT<double> fft;

   // detrend or demean
   x = detrend_and_demean_list(x, detrend, demean);

   for (size_t i = 0; i < n; ++i) {
     x_tmp = x[i];
     x_tmp.array() *= taper_array;
     x_padded = pad_vector(x_tmp, n_row, n_new);
     fft.fwd(x_fft, x_padded, 0);

     x_fft_list[i] = x_fft.array() * scale;
   }

   Eigen::VectorXcd t1;
   Eigen::VectorXcd t2;

   for (size_t i = 0; i < n; ++i) {
     t1 = x_fft_list[i];

     for (size_t j = i; j < n; ++j) {
       t2 = x_fft_list[j];
       t2 = t1.array() * t2.conjugate().array();

       // interpolate first value (do we want to do this?)
       t2[0] = 0.5 * (t2(1) + t2(n_new - 1));
       pgram[ind] = t2;
       ind += 1;
     }
   }


   // kernel multiplication
   if (spans(0) > 1) {
     Eigen::VectorXd kernel = modified_daniell(spans / 2);
     pgram = kernel_apply_list(pgram, kernel);
   }


   return(pgram);
 }
//==============================================================================



//******************************************************************************
// Welch
//******************************************************************************

//==============================================================================
//' @title
 //' spec_welch
 //'
 //' @description
 //' Calculate the periodogram using Welch's method.  This method only keeps the
 //' columns necessary for the transfer function calculation. This method is
 //' based on `spec.pgram`.
 //'
 //' @inheritParams spec.pgram
 //' @param length_subset length of each subset (integer)
 //' @param overlap percent to overlap subsets (double)
 //' @param window vector of length length_subset (numeric vector)
 //'
 //'
 //' @return periodogram from an input matrix using a Fast Fourier Transform and
 //' Welch's method.
 //'
 //' @noRd
 //'
 // [[Rcpp::export]]
 Eigen::MatrixXcd spec_welch(Eigen::MatrixXd& x,
                             size_t length_subset,
                             double overlap,
                             Eigen::VectorXd window
 ) {


   size_t n_row = x.rows();
   size_t n_col = x.cols();

   if (window.size() == 0) {
     VectorXd window = window_rectangle(length_subset);
   }
   if (window.size() != length_subset) {
     Rcpp::stop("spec_welch: `length_subset` must be the same length as window");
   }
   // check inputs
   if (length_subset > n_row) {
     Rcpp::stop("spec_welch: `length_subset` cannot be greater than the number of rows in `x`");
   }
   if (length_subset <= 0) {
     Rcpp::stop("spec_welch: `length_subset` must be a positive integer");
   }
   if (overlap > 1.0 || overlap < 0.0) {
     Rcpp::stop("spec_welch: `overlap` must be less than 1.0 and greater than 0.0");
   }


   size_t n_overlap = floor(length_subset * overlap);

   size_t n_new = next_n_eigen(length_subset);

   size_t n_fft = ceil((double)(n_row - length_subset) / (double)(length_subset - n_overlap));

   // suppose you want 10
   // size_t n_fft_10 = 1 - ceil((double)(n_row - length_subset)) / 10 / length_subset;

   // Rcpp::Rcout << "The value n_fft " << n_fft << std::endl;

   // get the starting index of each subset
   VectorXi starts = Eigen::VectorXi::LinSpaced(n_fft, 0, n_fft-1);
   starts = starts.array() * (length_subset - n_overlap);


   // window
   double scale = window_scale(window, n_new, n_fft);
   size_t n_col_pgram = n_col * double(n_col + 1) / 2;

   // temporary matrices and output
   MatrixXd x_block(n_col, length_subset);
   MatrixXcd x_fft_mat(n_new, n_col);
   MatrixXcd pgram_mat(n_new, n_col_pgram);
   pgram_mat.setZero();

   size_t ind;

   // calculate upper right -- doing in parallel requires additional work --
   for (size_t k = 0; k < n_fft; ++k) {
     ind = 0;
     x_block = x.middleRows(starts[k], length_subset);
     x_block = x_block.array().colwise() * window.array();
     x_block = detrend_and_demean_matrix(x_block, true, true);

     // calculate FFT
     x_fft_mat = fft_matrix(x_block, n_new);

     for (size_t i = 0; i < n_col; ++i) {
       for (size_t j = i; j < n_col; ++j) {

         pgram_mat.col(ind) = pgram_mat.col(ind).array() +
           x_fft_mat.col(i).array() * x_fft_mat.col(j).conjugate().array();


         ind += 1;

       }
     }
   }

   // pgram_mat = fill_lower_left(pgram_mat, n_col, start);

   pgram_mat.row(0) = 0.5 * (pgram_mat.row(1).array() + pgram_mat.row(n_new-1).array());
   pgram_mat *= scale;

   // MatrixXd m =  pgram_mat.cwiseAbs2();
   // double mx = m.col(2).maxCoeff();
   // Rcpp::Rcout << "The value max " << mx << std::endl;
   // std::complex<double> cplx_zero = std::complex<double>(1e-12, 1e-12);
   //
   // for (size_t i = 0; i < n_new; ++i) {
   //   if (abs(pgram_mat(i,2)) < 0.01) {
   //     pgram_mat(i, 2) = cplx_zero;
   //   }
   // }

   return(pgram_mat);
 }
//==============================================================================


// //==============================================================================
// //' @title
// //' spec_welch_trunc
// //'
// //' @description
// //' Calculate the periodogram using Welch's method.  This method only keeps the
// //' columns necessary for the transfer function calculation. This method is
// //' based on `spec.pgram`.
// //'
// //' @inheritParams spec_welch
// //'
// //'
// //' @return periodogram from an input matrix using a Fast Fourier Transform and
// //' Welch's method.
// //'
// //' @noRd
// //'
// // [[Rcpp::export]]
// Eigen::MatrixXcd spec_welch_trunc(const Eigen::MatrixXd& x,
//                                   size_t length_subset,
//                                   double overlap,
//                                   Eigen::VectorXd window
// ) {
//
//   return(spec_welch(x, length_subset, overlap, window, true));
//
// }
// //==============================================================================
//
//
// //==============================================================================
// //' @title
// //' spec_welch_complete
// //'
// //' @description
// //' Calculate the periodogram using Welch's method.  This method only keeps all
// //' columns even those that aren't necessary for the transfer function
// //' calculation. This method is based on `spec.pgram`.
// //'
// //' @inheritParams spec.pgram
// //'
// //'
// //' @return periodogram from an input matrix using a Fast Fourier Transform and
// //' Welch's method.
// //'
// //' @noRd
// //'
// // [[Rcpp::export]]
// Eigen::MatrixXcd spec_welch_complete(const Eigen::MatrixXd& x,
//                                      size_t length_subset,
//                                      double overlap,
//                                      Eigen::VectorXd window
// ) {
//
//   return(spec_welch(x, length_subset, overlap, window, false));
//
// }
// //==============================================================================

//******************************************************************************
// Estimate transfer functions
//******************************************************************************


//==============================================================================
//' @title
 //' solve_cplx_parallel
 //'
 //' @description
 //' Calculate the transfer function from a periodogram.
 //'
 //' @inheritParams spec.pgram
 //' @inheritParams make_groups
 //'
 //' @return the transfer functions.
 //'
 //' @noRd
 //'
 // [[Rcpp::export]]
 Eigen::MatrixXcd solve_cplx_parallel(const Eigen::MatrixXcd& x) {

   // original number of series
   size_t n_col = get_column_number(x.cols());

   size_t n_row = x.rows();
   size_t sub_size = n_col - 1;
   MatrixXcd out(n_row, sub_size);

   RcppThread::parallelFor(0, n_row, [&] (size_t i) {

     size_t ind = n_col;
     RowVectorXcd sub_v = x.row(i);

     VectorXcd y = sub_v.segment(1, sub_size).conjugate();

     MatrixXcd X(sub_size, sub_size);
     // const int p = sub_size;

     for (size_t i = 0; i < sub_size; ++i) {
       for (size_t j = i; j < sub_size; ++j) {

         X(i, j) = sub_v(ind);

         if (i != j) {
           X(j, i) = std::conj(sub_v(ind));
         }

         ind += 1;

       }
     }


     // ldlt/llt will be faster but sacrifices precision
     // Eigen::LLT<Eigen::MatrixXcd> llt(Eigen::MatrixXcd(p, p).setZero().selfadjointView<Lower>().
     //                                      rankUpdate(X.adjoint()));
     // out.row(i) = (llt.solve(X.adjoint() * y));


     // fullPivLu works fast for small matrices
     // ldlt/llt will be faster but sacrifices precision
     // https://eigen.tuxfamily.org/dox/group__TutorialLinearAlgebra.html
     out.row(i) = X.colPivHouseholderQr().solve(y);

   });

   return(out);
 }
//==============================================================================

// [[Rcpp::export]]
Eigen::MatrixXd list_to_matrix(Rcpp::List& x, int sub_start, int sub_end) {

  int n_cols = x.size();
  int n_rows = sub_end - sub_start + 1;

  Eigen::VectorXd sub_vector(n_rows);
  Eigen::VectorXd tmp;
  Eigen::MatrixXd y(n_rows, n_cols);


  for (int i = 0; i < n_cols; ++i) {
    tmp = x[i];
    y.col(i) = tmp.segment(sub_start, n_rows);
  }

  return(y);

}

//==============================================================================
//' @title
 //' solve_cplx_irr
 //'
 //' @description
 //' Calculate the transfer function from a periodogram with irregular sized
 //' groups. This is experimental to see if we can improve efficiency.
 //' Instead of fitting every frequency it fits groups of frequencies
 //' The goal is to lump many high frequency signals to increase signal to
 //' noise ratios, and only few low frequency signals to keep resolution at low
 //' frequency.
 //'
 //' @inheritParams make_groups
 //' @inheritParams fill_lower_left
 //'
 //'
 //' @return the transfer functions.
 //'
 //' @noRd
 //'
 // [[Rcpp::export]]
 Eigen::MatrixXcd solve_cplx_irr(Eigen::MatrixXcd& x,
                                 size_t n_groups) {

   // original number of series
   size_t n_col = get_column_number(x.cols());

   size_t n_row = x.rows() / 2 + 1; // half spectrum only
   VectorXi groups = make_groups(n_groups, n_row);
   size_t n_ols = groups.size();
   size_t sub_size = n_col - 1;

   VectorXi ind_sum(n_ols);

   ind_sum(0) = 0;
   for (size_t j = 1; j < n_ols; ++j) {
     ind_sum(j) = ind_sum(j - 1) + groups(j - 1);
   }

   MatrixXcd out(n_ols, sub_size);

   RcppThread::parallelFor(0, n_ols, [&] (size_t i) {

     size_t group_size = groups(i);

     MatrixXcd sub = x.middleRows(ind_sum(i), group_size);
     VectorXcd y(group_size * sub_size);

     MatrixXcd X(group_size * sub_size, sub_size);
     const int p = sub_size;

     // handle upper right
     size_t ind = n_col;

     for (size_t i = 0; i < sub_size; ++i) {
       for (size_t j = i; j < sub_size; ++j) {

         if (i == 0) {
           y.segment(j * group_size, group_size) = sub.col(ind - sub_size).conjugate();
         }

         X.col(j).segment(i * group_size, group_size) = sub.col(ind);

         if (i != j) {
           X.col(i).segment(j * group_size, group_size) = sub.col(ind).conjugate();
         }

         ind = ind + 1;
       }
     }

     // Eigen::LLT<Eigen::MatrixXcd> llt(Eigen::MatrixXcd(p, p).setZero().selfadjointView<Lower>().
     //                                      rankUpdate(X.adjoint()));
     // out.row(i) = (llt.solve(X.adjoint() * y));

     out.row(i) = X.colPivHouseholderQr().solve(y);
   });

   return(out);
 }





// //==============================================================================
// //' @title
// //' solve_cplx_irr2
// //'
// //' @description
// //' Calculate the transfer function from a periodogram with irregular sized
// //' groups. This is experimental to see if we can improve efficiency.
// //' Instead of fitting every frequency it fits groups of frequencies
// //' The goal is to lump many high frequency signals to increase signal to
// //' noise ratios, and only few low frequency signals to keep resolution at low
// //' frequency.
// //'
// //' @inheritParams make_groups
// //' @inheritParams fill_lower_left
// //'
// //'
// //' @return the transfer functions.
// //'
// //' @noRd
// //'
// // [[Rcpp::export]]
// Eigen::MatrixXcd solve_cplx_irr2(Eigen::MatrixXcd& x,
//                                 size_t n,
//                                 size_t max_lag) {
//
//   // original number of series
//   size_t n_col = get_column_number(x.cols());
//
//   size_t n_row = x.rows() / 2 + 1; // half spectrum only
//   ArrayXd groups = log_lags_eigen(n, max_lag-1);
//   size_t n_ols = groups.size();
//   size_t sub_size = n_col - 1;
//
//   VectorXi ind_sum(n_ols);
//
//   ind_sum(0) = 0;
//   for (size_t j = 1; j < n_ols; ++j) {
//     ind_sum(j) = ind_sum(j - 1) + groups(j - 1);
//   }
//
//   MatrixXcd out(n_ols, sub_size);
//
//   RcppThread::parallelFor (0, n_ols, [&] (size_t i) {
//
//     size_t group_size = groups(i);
//
//     MatrixXcd sub = x.middleRows(ind_sum(i), group_size);
//     VectorXcd y(group_size * sub_size);
//
//     MatrixXcd X(group_size * sub_size, sub_size);
//
//     // handle upper right
//     size_t ind = n_col;
//
//     for (size_t i = 0; i < sub_size; ++i) {
//       for (size_t j = i; j < sub_size; ++j) {
//
//         if (i == 0) {
//           y.segment(j * group_size, group_size) = sub.col(ind - sub_size).conjugate();
//         }
//
//         X.col(j).segment(i * group_size, group_size) = sub.col(ind);
//
//         if (i != j) {
//            X.col(i).segment(j * group_size, group_size) = sub.col(ind).conjugate();
//         }
//
//         ind = ind + 1;
//       }
//     }
//
//     out.row(i) = X.colPivHouseholderQr().solve(y);
//   });
//
//   return(out);
// }

//==============================================================================




//==============================================================================
//' @title
 //' ordinary_coherence_phase
 //'
 //' @description
 //' Calculate ordinary coherence and phase from a pgram. Reference:
 //' https://vru.vibrationresearch.com/lesson/coherence-mathematics/
 //'
 //' @param x periodogram matrix (complex matrix)
 //'
 //' @return Matrix with ordinary coherence and phase.
 //'
 //' @noRd
 //'
 // [[Rcpp::export]]
 Eigen::MatrixXd ordinary_coherence_phase(const Eigen::ArrayXXcd& x) {

   // original number of series
   size_t n_col = get_column_number(x.cols());

   if (n_col == 1) {
     Rcpp::stop("Cannot calculate coherency for a single periodogram");
   }

   size_t sub_size = n_col - 1;
   size_t n_coh_phase = n_col * (sub_size) / 2;
   size_t n_row = x.rows();
   size_t ind_0, ind_1, ind_2, ind_3;
   ind_0 = 0;

   MatrixXd coh_phase(n_row, n_coh_phase * 2);
   VectorXi diag(n_col);

   diag(0) = 0;
   for (size_t k = 0; k < sub_size; ++k) {
     diag(k+1) = diag(k) + n_col - k;
   }

   for (size_t i = 0; i < sub_size; ++i) {
     for (size_t j = i + 1; j < n_col; ++j) {

       ind_1 = diag(i) + j - i;
       ind_2 = diag(i);
       ind_3 = diag(j);

       coh_phase.col(ind_0) = (x.col(ind_1).abs2() /
         (x.col(ind_2) * x.col(ind_3))).real();

       coh_phase.col(ind_0 + n_coh_phase) = x.col(ind_1).arg();

       ind_0 += 1;
     }
   }

   // 0,1,2
   // 3,4,5
   // 6,7,8
   //
   // 0,1,2
   // -,3,4
   // -,-,5
   //
   // 0,1,2,3
   // -,4,5,6
   // -,-,7,8
   // -,-,-,9

   return(coh_phase);
 }
//==============================================================================



//******************************************************************************
// Convert from downsampled transfer functions to impulse response functions
// - Interpolate the downsampled transfer function to a smoothed response
// - apply the inverse FFT
//******************************************************************************




// //==============================================================================
// //' @title
// //' frequency_to_time_domain
// //'
// //' @description
// //' Convert from frequency to time domain
// //'
// //' @param pgram. (complex matrix)
// //' @param n_groups the frequency response functions. (numeric vector)
// //' @param knots knot positions for cubic spline interpolation. (numeric vector)
// //'
// //'
// //' @return the time domain cumulative impulse response.
// //'
// //' @noRd
// //'
// // [[Rcpp::export]]
// Eigen::MatrixXd frequency_to_time_domain(Eigen::MatrixXcd& pgram,
//                                          size_t n_groups
//                                          ) {
//
//   // Get frequencies
//   ArrayXd fs_all = determine_frequency(pgram.rows());
//   ArrayXd fs_sub = group_frequency(fs_all, n_groups);
//
//   // Determine transfer function
//   MatrixXcd tf = solve_cplx_irr(pgram, n_groups);
//
//   // Determine knot spacing
//   VectorXd knots = power_spaced(n_groups / 2, fs_all.minCoeff(), fs_all.maxCoeff(), 10);
//
//   // Smoothly interpolate transfer function
//   MatrixXcd tf_i = interpolate_tf(tf, fs_sub, fs_all, knots);
//
//   // Perform inverse fft on transfer function
//   Eigen::FFT<double> fft;
//   size_t n_row = tf_i.rows();
//   size_t n_half = n_row / 2;
//   size_t n_col = tf_i.cols();
//
//   VectorXcd tf_i_out(n_row);
//   VectorXd z(n_row);
//   MatrixXd out(n_half, n_col);
//
//   // remove DC components of transfer function
//   tf_i.row(0) = tf_i.row(1);
//   tf_i.row(n_half) = tf_i.row(n_half - 1);
//
//   for (size_t i = 0; i < n_col; ++i){
//     tf_i_out.setZero();
//     z.setZero();
//     // x_out.head(n_row) = x.col(i);
//     tf_i_out = tf_i.col(i);
//     fft.inv(z, tf_i_out);
//     out.col(i) = z.head(n_half).array() + z.reverse().head(n_half).array();
//   }
//
//   return(out);
// }
// //==============================================================================


// frequency_to_time_domain <- function(x, y) {
//
//   dc1 <- y[1,] # DC values
//   dc2 <- y[nrow(y),] # DC values
//   x_n <- x[-c(1, length(x))] # DC indices
//   y_n <- y[-c(1, nrow(y)),] # DC values
//   knots <- hydrorecipes:::log_lags(20, max(x_n))
//   knots <- knots[-c(1, length(knots))]
//   len <- min(knots):max(knots)
//
// # this is used to smooth the complex response
//   sp_in  <- splines2::bSpline(x_n, knots = knots, Boundary.knots = c(0, max(x_n)+1e-16))
//     sp_out <- splines2::bSpline(len, knots = knots, Boundary.knots = c(0, max(x_n)+1e-16))
//
//
//     out <- matrix(NA_real_,
//                   nrow = length(len),
//                   ncol = ncol(y_n) + 1)
//     out[, 1] <- len
//
// # loop through each transfer function - smooth the response - estimate brf
//     for (i in 1:ncol(y_n)) {
//
// # smooth the real and imaginary components separately
//       re <- Re(y_n[, i])
//         im <- Im(y_n[, i])
//
//         fit_r <- RcppEigen::fastLm(X = sp_in, y = re)$coefficients
//         fit_i <- RcppEigen::fastLm(X = sp_in, y = im)$coefficients
//
//         re <- as.vector(sp_out %*% fit_r)
//         im <- as.vector(sp_out %*% fit_i)
//
// # estimate the brf from the smoothed uniform spaced frf
//         out[, i + 1] <- brf_from_frf(complex(real = re, imaginary = im),
//                                      dc1[1], dc2[1])
//     }
//
//
//     out
// }


// //==============================================================================
// //' @title
// //' frf_to_brf
// //'
// //' @description
// //' Convert from frequency to time domain.  This function uses the half
// //' frequency response.
// //'
// //' @param x the frequency response function. (complex matrix)
// //' @param dc1 the first value of the FRF. (double)
// //' @param dc2 the last value of the half FRF. (double)
// //'
// //'
// //'
// //' @return the time domain cumulative impulse response.
// //'
// //' @noRd
// //'
// // [[Rcpp::export]]
// Eigen::ArrayXd frf_to_brf(const Eigen::VectorXcd& x,
//                           std::complex<double> dc1,
//                           std::complex<double> dc2) {
//
//
//   Eigen::FFT<double> fft;
//   size_t n = x.size();
//   Eigen::VectorXcd x_out(n * 2 + 2);
//   x_out << dc1, x, dc2, x.reverse().conjugate();
//
//   size_t n_full = x_out.size();
//   Eigen::VectorXcd z(n_full);
//   Eigen::ArrayXd out(n_full);
//
//   fft.inv(z, x_out);
//   out = z.real();
//
//   out = (out.head(n) + out.reverse().head(n));
//
//   return(out);
//
// }
// //==============================================================================


// //==============================================================================
// // x are frequencies
// // y is the transfer function
// // knots for the spline regression fit
// // degree fro the spline regresssion fit
// //' @title
// //' interpolate_tf
// //'
// //' @description
// //' Go from irregularly spaced frequency response function to a regularly
// //' spaced version.  This function will smooth out local variability depending
// //' on the chosen knots.
// //'
// //' @inheritParams which_indices
// //'
// //' @param x the frequencies. (numeric vector)
// //' @param y the frequency response function values. (complex matrix)
// //' @param knots locations used to interpolate the frequency response. (numeric vector)
// //' @param degree the degree for the \code{b_spline} function. (integer)
// //' @param x_interp the frequencies for interpolation. (numeric vector)
// //'
// //'
// //'
// //' @return the time domain cumulative impulse response.
// //'
// //' @noRd
// //'
// // [[Rcpp::export]]
// Eigen::MatrixXcd interpolate_tf(Eigen::MatrixXcd& x,
//                                  const Eigen::ArrayXd& frequency_irregular,
//                                  const Eigen::ArrayXd& frequency_regular,
//                                  Eigen::VectorXd& knots
// ) {
//
//   size_t n_x = x.rows();   // transfer function values
//   size_t n_freq = frequency_irregular.size(); // frequencies for transfer function
//   size_t n_col = x.cols(); // number of series
//   size_t n = n_x - 2;      // n without dc components
//   size_t n_knots = knots.size();
//
//
//   // check inputs
//   if (n_x != n_freq) stop("interpolate_tf: x and y lengths must be equal");
//   if (n_x < n_knots) stop("interpolate_tf: the number of knots cannot be greater than the input length");
//
//
//   // check knots
//   for (size_t i = 0; i < knots.size(); ++i) {
//     if (knots(i) < frequency_irregular(i + 1)) {
//       knots(i) = (frequency_irregular(i + 1) + frequency_irregular(i + 2)) / 2.0;
//     }
//   }
//
//   // we don't want to use the dc values in the interpolation
//   VectorXcd dc1 = x.row(0);
//   VectorXcd dc2 = x.row(n_x - 1);
//
//
//   VectorXd frequency_sub = frequency_irregular.segment(1, n);
//   size_t max = frequency_regular.size();
//
//   // generate knots
//   MatrixXd bs_in  = b_spline(frequency_sub, knots);
//   MatrixXd bs_out = b_spline(frequency_regular, knots);
//   // size_t n_interp = x_interp.size();
//
//   VectorXcd frf(max);
//   VectorXd fit_real(n);
//   VectorXd fit_imag(n);
//   VectorXd re(n);
//   VectorXd im(n);
//   VectorXcd x_sub(n);
//   Eigen::MatrixXcd out(max * 2 , n_col);
//
//
//   for (size_t i = 0; i < n_col; ++i) {
//
//     // do y by col
//     x_sub = x.col(i).segment(1, n);
//
//     fit_real = bs_in.colPivHouseholderQr().solve(x_sub.real());;
//     fit_imag = bs_in.colPivHouseholderQr().solve(x_sub.imag());;
//
//     re = bs_out * fit_real;
//     im = bs_out * fit_imag;
//
//     // generate full length sequence
//     for (size_t j = 0; j < max; ++j) {
//       frf(j) = std::complex<double>(re[j], im[j]);
//     }
//       // Rcpp::Rcout << "The value frf head " << frf.head(5) << std::endl;
//       // Rcpp::Rcout << "The value frf tail" << frf.head(5) << std::endl;
//
//     // out.col(i) << dc1(i), frf, dc2(i), frf.reverse().conjugate();
//     // out.col(i) << y_sub(1), frf, y_sub.tail(1), frf.reverse().conjugate();
//     out.col(i) << frf, frf(0), frf.tail(max-1).reverse().conjugate();
//     // out.col(i) << frf(0), frf, frf.tail(1), frf.reverse().conjugate();
//       // Rcpp::Rcout << "The value out head " << out.col(i).head(5) << std::endl;
//       // Rcpp::Rcout << "The value out tail" << out.col(i).tail(5) << std::endl;
//
//
//   }
//
//   return(out);
//
// }
// //==============================================================================


//==============================================================================
//' @title
 //' transfer_pgram_smooth
 //'
 //' @description
 //' Calculate the transfer function from an input matrix. This function uses
 //' irregular sized groups using `make_groups`. This is experimental to see if
 //' and designed to be relatively fast. Instead of fitting every frequency and
 //' aggregating post solving, it fits groups of frequencies.
 //' The goal is to lump many high frequency signals to increase signal to
 //' noise ratios, and only few low frequency signals to keep resolution at low
 //' frequency.
 //'
 //' @inheritParams spec.pgram
 //' @inheritParams make_groups
 //' @param n_col number of covariate columns (integer)
 //'
 //'
 //' @return the transfer functions.
 //'
 //' @noRd
 //'
 // [[Rcpp::export]]
 Eigen::MatrixXcd transfer_pgram_smooth(Eigen::MatrixXd& x,
                                        const Eigen::VectorXi& spans,
                                        bool detrend,
                                        bool demean,
                                        double taper,
                                        double power,
                                        size_t n_groups) {

   // size_t min_aggregate = 1;

   // if (spans(0) == 0) {
   //   min_aggregate = 3;
   // }

   MatrixXcd pgram = spec_pgram(x, spans, detrend, demean, taper);
   MatrixXcd out   = solve_cplx_irr(pgram, n_groups);

   return(out);
 }
//==============================================================================


//==============================================================================
// [[Rcpp::export]]
Eigen::MatrixXcd transfer_pgram(Eigen::MatrixXd& x,
                                const Eigen::VectorXi& spans,
                                bool detrend,
                                bool demean,
                                double taper) {

  Eigen::MatrixXcd pgram = spec_pgram(x, spans, detrend, demean, taper);
  Eigen::MatrixXcd out   = solve_cplx_parallel(pgram);
  return(out);
}
//==============================================================================




//==============================================================================
//' @title
 //' transfer_welch
 //'
 //' @description
 //' Calculate the transfer function from an input matrix. This function uses
 //' irregular sized groups using `make_groups`. This is experimental to see if
 //' and designed to be relatively fast. Instead of fitting every frequency and
 //' aggregating post solving, it fits groups of frequencies.
 //' The goal is to lump many high frequency signals to increase signal to
 //' noise ratios, and only few low frequency signals to keep resolution at low
 //' frequency.
 //'
 //' @inheritParams spec.pgram
 //' @inheritParams make_groups
 //' @param n_col number of covariate columns (integer)
 //'
 //'
 //' @return the transfer functions.
 //'
 //' @noRd
 //'
 // [[Rcpp::export]]
 Eigen::MatrixXcd transfer_welch(Eigen::MatrixXd& x,
                                 size_t length_subset,
                                 double overlap,
                                 Eigen::VectorXd window
 ) {

   MatrixXcd pgram = spec_welch(x, length_subset, overlap, window);
   MatrixXcd out   = solve_cplx_parallel(pgram);

   return(out);
 }
//==============================================================================


// //==============================================================================
// //' @title
// //' predict_pgram_frf
// //'
// //' @description
// //' Calculate the transfer function from an input matrix. This function uses
// //' irregular sized groups using `make_groups`. This is experimental to see if
// //' and designed to be relatively fast. Instead of fitting every frequency and
// //' aggregating post solving, it fits groups of frequencies.
// //' The goal is to lump many high frequency signals to increase signal to
// //' noise ratios, and only few low frequency signals to keep resolution at low
// //' frequency.
// //'
// //' @inheritParams spec.pgram
// //' @inheritParams make_groups
// //' @param n_col number of covariate columns (integer)
// //'
// //'
// //' @return the transfer functions.
// //'
// //' @noRd
// //'
// // [[Rcpp::export]]
// Eigen::MatrixXd predict_pgram_frf(Eigen::MatrixXd& x,
//                                   Eigen::MatrixXd& x_out,
//                                   Eigen::VectorXi spans, // spec_pgram
//                                   size_t n_groups) {
//
//   size_t n = x.rows();
//   size_t n_out = x_out.rows();
//
//   // Get frequencies
//   ArrayXd fs_all = determine_frequency(n);
//   ArrayXd fs_out = determine_frequency(n_out);
//   ArrayXd fs_sub = group_frequency(fs_all, n_groups);
//
//   // calculate periodogram
//   MatrixXcd pgram = spec_pgram(x, spans, true, true, 0.1);
//
//   // estimate transfer function
//   MatrixXcd tf = solve_cplx_irr(pgram, n_groups);
//
//   // Determine knot spacing
//   VectorXd knots = power_spaced(n_groups / 2, fs_all.minCoeff(), fs_all.maxCoeff(), 10);
//
//   // interpolate to length similar to x_out
//   MatrixXcd x_interp = interpolate_tf(tf, fs_sub, fs_out, knots);
//
//   MatrixXd ret(x_interp.rows(), x_interp.cols());
//
//   Rcpp::Rcout << "The value frf rows " << tf.rows() << std::endl;
//   Rcpp::Rcout << "The value pgram rows " << pgram.rows() << std::endl;
//   Rcpp::Rcout << "The value x_interp rows " << x_interp.rows() << std::endl;
//   Rcpp::Rcout << "The value x_out rows " << x_out.rows() << std::endl;
//
//   for (size_t i = 0; i < x_interp.cols(); ++i) {
//     ret.col(i) = convolve_tf(x_out.col(i+1),
//                              x_interp.col(i).reverse());
//   }
//
//   return(ret);
// }
// //==============================================================================



//
// //==============================================================================
// // [[Rcpp::export]]
// Eigen::MatrixXd predict_pgram_smooth_frf(Eigen::MatrixXd& x,
//                                          Eigen::VectorXi span,
//                                          double power,
//                                          size_t n_groups) {
//
//
//   size_t n_col = x.cols();
//   size_t min_aggregate = 1;
//   Rcpp::Rcout << "The value x " << x.rows() << std::endl;
//   Eigen::MatrixXcd pgram = spec_pgram(x, span, false, false, 0.1);
//   Rcpp::Rcout << "The value pgram " << pgram.rows() << std::endl;
//   Rcpp::Rcout << "The value pgram cols " << pgram.cols() << std::endl;
//
//   size_t n_row = pgram.rows() / 2 + 1; // half spectrum only
//
//   VectorXd fs_all = determine_frequency(n_row - 2);
//   VectorXd fs_sub = group_frequency(fs_all, power, n_groups, min_aggregate);
//   // double fs_diff = fs_all(1) - fs_all(0);
//   Eigen::ArrayXd knots =power_spaced(20, fs_all.minCoeff(), fs_all.maxCoeff(), 2.73);
//   // knots(0) = fs_all.minCoeff();
//   // knots(knots.size()-1) = fs_all.maxCoeff();
//   // Eigen::ArrayXd knots(11);
//   // knots << fs_all.minCoeff(),
//   //          1.157407e-05/2, 1.157407e-05, 1.157407e-05*1.5, 1.157407e-05*2,
//   //          1.157407e-05*4, 1.157407e-05*10, 1.157407e-05*100, 1.157407e-05*1000,
//   //          1.157407e-05*10000,
//   //          fs_all.maxCoeff();
//   Eigen::MatrixXcd out = solve_cplx_irr(pgram,
//                                               power,
//                                               n_groups,
//                                               min_aggregate);
//
//   // for (size_t i = 0; i < fs_sub.size()/4; ++i){
//   //   knots(i) = fs_sub.segment(i*4, 3).mean();
//   // }
//   // knots(0) = fs_all(0);
//   // knots(fs_sub.size()/4) = fs_all(pgram.rows() / 2 - 2);
//   Eigen::MatrixXcd frf = interpolate_tf(fs_sub, out, knots, 3, fs_all);
//
//   Eigen::MatrixXd ret(frf.rows(), frf.cols());
//   // size_t sub_size = frf.cols();
//
//   Rcpp::Rcout << "The value frf rows " << frf.rows() << std::endl;
//   Rcpp::Rcout << "The value pgram rows " << pgram.rows() << std::endl;
//   Rcpp::Rcout << "The value x rows " << x.rows() << std::endl;
//
//   for (size_t i = 1; i < frf.cols(); ++i) {
//     size_t tt = 5;
//     for (size_t j = 0; j < pgram.rows(); ++j) {
//       if (std::abs(pgram(j, tt)) < 1e4) {
//         frf(j, i) = std::complex<double>(0.0, 0.0);
//       }
//     }
//   }
//
//   for (size_t i = 0; i < frf.cols(); ++i) {
//     ret.col(i) = convolve_tf(
//       pad_vector(
//         x.col(i+1),
//         x.rows(),
//         frf.rows()),
//         frf.col(i).reverse());
//   }
//
//
//   // Eigen::FFT<double> fft;
//   // size_t n = irf.size();
//   // Eigen::VectorXcd x_out(n * 2 + 2);
//   // Eigen::VectorXcd dc1 = out.row(0);
//   // Eigen::VectorXcd dc2 = out.row(out.rows() - 1);
//   //
//   // // x_out << dc1, x, dc2, x.reverse().conjugate();
//   // Rcpp::Rcout << "The value xout " << x_out.rows() << std::endl;
//
//
//   return(ret);
// }




// //==============================================================================
// //' @title
// //' solve_weighted_cplx_irr
// //'
// //' @description
// //' Calculate the transfer function from a periodogram with irregular sized
// //' groups. This is experimental to see if we can improve efficiency.
// //' Instead of fitting every frequency it fits groups of frequencies
// //' The goal is to lump many high frequency signals to increase signal to
// //' noise ratios, and only few low frequency signals to keep resolution at low
// //' frequency. This function weights the central frequency the most in the
// //' group.
// //'
// //' @inheritParams spec.pgram
// //' @inheritParams make_groups
// //' @param n_col number of covariate columns (integer)
// //'
// //'
// //' @return the transfer functions.
// //'
// //' @noRd
// //'
// // [[Rcpp::export]]
// Eigen::MatrixXcd solve_weighted_cplx_irr(Eigen::MatrixXcd& x,
//                                                Eigen::VectorXi groups,
//                                                Eigen::VectorXi sizes,
//                                                size_t n_col) {
//
//   size_t sub_size = n_col - 1;
//
//   size_t n_groups = groups.size();
//
//   MatrixXcd out(n_groups, sub_size);
//
//   RcppThread::parallelFor (0, n_groups, [&] (size_t i) {
//     size_t group_size = sizes(i);
//     MatrixXcd sub = x.middleRows(groups(i), group_size);
//     VectorXcd y(group_size * sub_size);
//     MatrixXcd X(group_size * sub_size, sub_size);
//     size_t rem;
//     size_t int_div;
//
//     for (size_t j = 0; j < n_col * sub_size; ++j) {
//       rem = j % n_col;
//       int_div = floor(j / n_col);
//
//       if (rem == 0) {
//         y.segment(int_div * group_size, group_size) = sub.col(j).array() * window_hann_eigen(group_size).array();
//       } else {
//         X.col(rem - 1).segment(int_div * group_size, group_size) = sub.col(j).array()* window_hann_eigen(group_size).array();
//       }
//     }
//
//     out.row(i) = X.householderQr().solve(y);
//   });
//
//
//   return(out);
// }
// //==============================================================================

// // [[Rcpp::export]]
// Eigen::VectorXd make_groups_test(size_t n_row,
//                             size_t power,
//                             size_t n_groups) {
//
//
//   Eigen::VectorXd groups_f = VectorXd::LinSpaced(n_groups + 2, 0, 1);
//   groups_f = groups_f.array().pow(power);
//   groups_f(0) = 0.0;
//   groups_f(n_groups + 1) = 0.0;
//   double scale = std::floor(n_row / groups_f.sum());
//
//   groups_f = ((groups_f * scale).array().max(1.0).round());
//   Eigen::VectorXi groups(n_groups+2);
//
//   for (size_t i = 0; i < n_groups+2; ++i) {
//     groups(i) = (size_t)groups_f(i);
//   }
//   return(groups);
// }



//
// // [[Rcpp::export]]
// Eigen::VectorXi make_groups_sigmoid(size_t n_row,
//                                     double min_sigmoid,
//                                     size_t n_groups) {
//
//   VectorXf groups_f = -VectorXf::LinSpaced(n_groups, min_sigmoid, 6);
//   groups_f = 1.0 / (1.0 + groups_f.array().exp());
//   double scale = std::floor(n_row / 2.0 / groups_f.sum());
//   groups_f = ((groups_f * scale).array().max(1).round());
//
//   Eigen::VectorXi groups(n_groups);
//   for (size_t i = 0; i < n_groups; ++i) {
//     groups(i) = (size_t)groups_f(i);
//   }
//   size_t g_sum = groups.sum();
//   if (g_sum > n_row/2) {
//     groups(n_groups - 1) -= (g_sum - n_row / 2);
//   }
//   return(groups);
//
// }

// // [[Rcpp::export]]
// Eigen::MatrixXcd so(Eigen::MatrixXcd& X,
//                             Eigen::VectorXcd& y) {
//
//   return(X.householderQr().solve(y));
// }
//==============================================================================


// // [[Rcpp::export]]
// Eigen::VectorXd basis_eigen(const Eigen::ArrayXd& x,
//                             size_t degree,
//                             size_t i,
//                             const Eigen::VectorXd& knots
// ) {
//
//   size_t n = x.size();
//   Eigen::VectorXd out = VectorXd::Zero(n);
//
//   Eigen::ArrayXd alpha1 = VectorXd::Zero(n);
//   Eigen::ArrayXd alpha2 = VectorXd::Zero(n);
//
//   if (degree == 0) {
//     out.setOnes();
//     //   for (size_t j = 0; j < n; j++) {
//     //     if (x(j) >= knots(i) & x(j) < knots(i+1)) {
//     //       out(j) = 1.0;
//     //     } else {
//     //       out(j) = 0.0;
//     //     }
//     //   }
//   } else {
//     if (knots(degree + i) - knots(i) != 0) {
//       alpha1 = (x.array() - knots(i)) / (knots(degree + i) - knots(i));
//     }
//     if (knots(degree + i + 1) - knots(i + 1) != 0) {
//       alpha2 = (knots(i + degree + 1) - x.array()) /
//         (knots(degree + i + 1) - knots(i + 1));
//     }
//     out = alpha1 * basis_eigen(x, (degree-1), i, knots).array() +
//       alpha2 * basis_eigen(x, (degree-1), (i+1), knots).array();
//   }
//
//   return(out);
// }
//
//
// // [[Rcpp::export]]
// Eigen::MatrixXd bs_eigen(const Eigen::ArrayXd& x,
//                          size_t degree,
//                          const Eigen::VectorXd& knots) {
//
//   size_t n_knots = knots.size();
//   size_t n_x = x.size();
//   size_t k = n_knots + degree - 1;
//
//   Eigen::MatrixXd out = MatrixXd::Zero(n_x, k);
//   Eigen::ArrayXd knots_pad(degree * 2 + n_knots);
//   Eigen::VectorXd pad(degree);
//
//   pad.setConstant(knots(0));
//   knots_pad.head(degree) = pad;
//   pad.setConstant(knots(n_knots-1));
//   knots_pad.tail(degree) = pad;
//   knots_pad.segment(degree, n_knots) = knots;
//
//
//   // for (size_t i = 0; i < k; ++i) {
//   RcppThread::parallelFor (0, k, [&] (size_t i) {
//
//     out.col(i) = basis_eigen(x, degree, i, knots_pad);
//     // }
//   });
//
//   for (size_t j = n_x-1; j > 0; --j) {
//     if (x(j) >= knots(n_knots-1)){
//       out(j, k-1) = 1.0;
//     } else {
//       break;
//     }
//   }
//
//   return(out);
// }
//
// // [[Rcpp::export]]
// Eigen::MatrixXi which_indices2(const Eigen::VectorXd& x,
//                                const Eigen::VectorXd& knots
// ) {
//
//   size_t n_x = x.size();
//   size_t n_knots = knots.size();
//   Eigen::MatrixXi ind(n_knots, 3);
//
//   double current_val;
//   double current_knot = knots(0);
//
//   size_t s = 0;
//   size_t n = 1;
//
//   for (size_t i = 0; i < n_knots; ++i) {
//     ind(i, 1) = (x.array() <= knots(i)).count();
//     ind(i, 0) = knots(i);
//   }
//   ind(n_knots, 2) = 0.0;
//   for (size_t i = 1; i < n_knots; ++i) {
//     ind(i-1, 2) = ind(i, 1) - ind(i-1, 1);
//   }
//
//
//   return(ind);
// }



// // [[Rcpp::export]]
// Eigen::VectorXcd spline_fit(Eigen::RowVectorXd x,
//                             Eigen::RowVectorXcd& points,
//                             Eigen::VectorXd& x_out) {
//
//   double min_x = x.row(0).array().minCoeff();
//   double max_x = x.row(0).array().maxCoeff();
//   size_t n_p = points.size();
//
//   MatrixXd points_re(2, n_p);
//   MatrixXd points_im(2, n_p);
//   // VectorXd knots = Eigen::VectorXd::LinSpaced(3,0,1);
//   points_re.row(0) = x;
//   points_im.row(0) = x;
//   points_re.row(1) = points.real();
//   points_im.row(1) = points.imag();
//
//   Eigen::Spline2d spline_re = Eigen::SplineFitting<Spline2d>::Interpolate(points_re, 1);
//   Eigen::Spline2d spline_im = Eigen::SplineFitting<Spline2d>::Interpolate(points_im, 1);
//
//   size_t n = x_out.size();
//   x_out = (x_out.array() - min_x) / (max_x - min_x);
//   Eigen::VectorXcd values(n);
//
//   // for (int i = 0; i < n; ++i) {
//   RcppThread::parallelFor (0, n, [&] (size_t i) {
//     double re = spline_re(x_out(i))[1];
//     double im = spline_im(x_out(i))[1];
//     std::complex<double> cplx(re, im);
//     values(i) = cplx;
//   });
//
//   return(values);
//
// }
//
// // [[Rcpp::export]]
// Eigen::RowVectorXd spline_fit2(Eigen::MatrixXd ctrls,
//                                Eigen::RowVectorXd knots_in,
//                                Eigen::VectorXd x_out) {
//
//   size_t n_col = ctrls.cols();
//   size_t n_knots = knots_in.size();
//   // size_t knot_buffer = (n_col - n_knots) / 2 + 1;
//   double min_x = ctrls.row(0).array().minCoeff();
//   double max_x = ctrls.row(0).array().maxCoeff();
//   // RowVectorXd knots = Eigen::RowVectorXd::Zero(n_knots + 2 * knot_buffer);
//   //
//   // RowVectorXd buffer_values(knot_buffer);
//   //
//   // knots.segment(knot_buffer, n_knots) = knots_in;
//   // knots = (knots.array() - min_x) / (max_x - min_x);
//
//   // ctrls.row(0) = (ctrls.row(0).array() - min_x) / (max_x - min_x);
//
//   // buffer_values.setConstant(1.0);
//   // knots.tail(knot_buffer) = buffer_values;
//   // buffer_values.setConstant(0.0);
//   // knots.head(knot_buffer) = buffer_values;
//
//   Eigen::Spline2d spline = Eigen::Spline<double, 2, Dynamic>(knots_in, ctrls);
//
//   size_t n = x_out.size();
//   x_out = (x_out.array() - min_x) / (max_x - min_x);
//   Eigen::MatrixXd values(n,2);
//   RowVectorXd re = spline.BasisFunctions(0.5, 3, knots_in);
//
//   // re = spline(0.351659507062997);
//   // for (int i = 0; i < n; ++i) {
//   // // RcppThread::parallelFor (0, n, [&] (size_t i) {
//   //   re = spline(x_out(i));
//   //   // std::complex<double> cplx(re, im);
//   //   values.row(i) = re;
//   // // });
//   // }
//
//   return(re);
// }
//
// typedef Eigen::Spline<double, 1, 2> Spline1D;
// typedef Eigen::SplineFitting<Spline1D> SplineFitting1D;
//
// // [[Rcpp::export]]
// Eigen::VectorXd spline_fit3(Eigen::RowVectorXd x,
//                             Eigen::RowVectorXd y,
//                             Eigen::VectorXd x_out) {
//
//
//   double min_x = x.array().minCoeff();
//   double max_x = x.array().maxCoeff();
//   size_t n = x_out.size();
//   x_out = (x_out.array() - min_x) / (max_x - min_x);
//   x = (x.array() - min_x) / (max_x - min_x);
//
//   const auto fit = SplineFitting1D::Interpolate(y, 2, x);
//   Spline1D spline(fit);
//   Eigen::VectorXd values(n);
//   double x_in;
//
//   for (size_t i = 0; i < n; ++i) {
//     x_in = x_out(i);
//     values(i) = spline(x_in).coeff(0);
//   }
//
//
//
//   return(values);
//
// }
//
//
// // [[Rcpp::export]]
// Eigen::VectorXd spline_fit4(Eigen::MatrixXd ctrls,
//                             Eigen::RowVectorXd knots_in,
//                             Eigen::VectorXd x_out) {
//
//   size_t n_col = ctrls.cols();
//   size_t n_knots = knots_in.size();
//   size_t knot_buffer = (n_col - n_knots) / 2 + 1;
//   double min_x = ctrls.row(0).minCoeff();
//   double max_x = ctrls.row(0).maxCoeff();
//   RowVectorXd knots = Eigen::RowVectorXd::Zero(n_knots + 2 * knot_buffer);
//
//   RowVectorXd buffer_values(knot_buffer);
//   knots_in = (knots_in.array() - min_x) / (max_x - min_x);
//   knots.segment(knot_buffer, n_knots) = knots_in;
//   // knots = (knots.array() - min_x) / (max_x - min_x);
//
//   ctrls.row(0) = (ctrls.row(0).array() - min_x) / (max_x - min_x);
//
//   buffer_values.setConstant(1.0);
//   knots.tail(knot_buffer) = buffer_values;
//   buffer_values.setConstant(0.0);
//   knots.head(knot_buffer) = buffer_values;
//
//   Eigen::Spline2d spline = Eigen::Spline<double, 2, Dynamic>(knots, ctrls);
//   Eigen::VectorXd test = spline.BasisFunctions(0.5, 3, knots_in);
//
//   return(test);
// }


// // [[Rcpp::export]]
// Eigen::VectorXd solve_test(Eigen::MatrixXd X,Eigen::VectorXd  y) {
//
//   return(X.householderQr().solve(y));
//
// }



//******************************************************************************


// //' @title
// //' fft_one
// //'
// //' @description
// //' Do an FFT on a vector
// //'
// //' @inheritParams fft_matrix
// //'
// //'
// //' @return the FFT of a vector with n_new values
// //'
// //' @noRd
// //'
// // [[Rcpp::export]]
// Eigen::VectorXcd fft_one(const Eigen::VectorXd& x,
//                          size_t n_new){
//
//   size_t n_row = x.rows();
//
//   Eigen::VectorXcd x_fft(n_new);
//   Eigen::VectorXd  x_padded = pad_vector(x, n_row, n_new);
//
//   Eigen::FFT<double> fft;
//
//   fft.fwd(x_fft, x_padded, 0);
//
//
//   return(x_fft);
// }
//
//

/*** R

x <- rnorm(1e7)
y <- rnorm(1e6)

# plot(frecipes:::convolve_overlap_save(x, rev(y)), type = 'l')
# plot(frecipes:::convolve_overlap_add(x, rev(y)), type = 'l', col = "green")
# points(frecipes:::convolve_filter(x, y, TRUE, TRUE), type = 'l', col = 'red')
y1 <- rev(y)
tmp <-bench::press(
  y_len = c(1e2, 3e5),
  {
    x <- rnorm(1e7)
    y <- rnorm(y_len)
    y1 <- rev(y)
    bench::mark(
      frecipes:::convolve_filter(x, y, TRUE, TRUE),
      frecipes:::convolve_overlap_add(x, y1),
      frecipes:::convolve_overlap_save(x, y1, 0),
      # frecipes:::convolve_overlap_save(x, y1, 1),
      # frecipes:::convolve_overlap_save(x, y1, 2),
      # frecipes:::convolve_vec(x, frecipes:::pad_vector(y1, y_len, 1e7)),
      check = TRUE,
      min_iterations = 2
    )}
)

tmp
x <- cumsum(rnorm(1000000))
y <- frecipes:::window_nuttall(1000)/sum(window_nuttall(1000))
plot(x, type = 'l')
points(frecipes:::convolve_overlap_save(x, y, 1), type = 'l', col = 'red')


m <- matrix(rep(y, 10), ncol = 10)
l <- list(y,y,y,y,y,y,
          # y,y,y,y,y,y,
          # y,y,y,y,y,y,
          y,y,y,y)
bench::mark(
  a <- frecipes:::convolve_overlap_save_list(x, l),
  b <- frecipes:::convolve_list(x, l, TRUE, TRUE),
  # d <- frecipes:::convolve_filter(x, y, TRUE, TRUE),
  # e <- frecipes:::convolve_matrix(x, m, TRUE, TRUE),
  check = FALSE,
  min_iterations = 1
)

bench::mark(
  a <- frecipes:::spec_pgram(m, spans = 3, TRUE, TRUE, taper = 0.1),
  base <- spec.pgram(m, spans = 3, demean = TRUE, detrend = TRUE, taper = 0.1, plot = FALSE),
  b <- collapse::qM(frecipes:::spec_pgram_list(l, spans = 3, TRUE, TRUE, taper = 0.1)),
  check = FALSE
)



*/
