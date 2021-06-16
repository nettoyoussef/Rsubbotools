
#include<RcppGSL.h>
// [[Rcpp::depends(RcppGSL)]]

#include "common.h"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>

//' Returns the Fisher Information matrix and its inverse for the AEP
//'
//' Returns the Fisher Information matrix and its inverse for the Asymmetric
//' Power Exponential distribution for the given parameters.
//'
//' @param size (numeric) - number of observations (Default: 01)
//' @param bl (numeric) - set the left exponent (Default: 2.0)
//' @param br (numeric) - set the right exponent (Default: 2.0)
//' @param m  (numeric) - the location parameter (Default: 0.0)
//' @param al (numeric) - the left scale parameter (Default: 1.0)
//' @param ar (numeric) - the right scale parameter (Default: 1.0)
//' @param O_munknown (numeric) - if true assumes \eqn{m} is known
//' @return a list containing three elements:
//' * std_error - the standard error for the parameters
//' * infmatrix - the Fisher Information Matrix
//' * inv_infmatrix - the Inverse Fisher Information Matrix
//' @export
//' @md
// [[Rcpp::export]]
Rcpp::List subboafish(
                       unsigned size = 1
                      ,double bl=2.0
                      ,double br=2.0
                      ,double m=0.0
                      ,double al=1.0
                      ,double ar=1.0
                      ,unsigned O_munknown=0
                      ){

  // dimension of the matrix
  size_t dim=(O_munknown?4:5);

  // name of parameters
  Rcpp::CharacterVector par(dim);

  // information matrix
  // allocates space
  RcppGSL::Matrix I(dim, dim);
  RcppGSL::Matrix invI(dim, dim);

  // needed for inversion
  RcppGSL::Matrix J(dim, dim);
  gsl_permutation *P;
  P = gsl_permutation_alloc (dim);
  int signum;

  // the order is bl,br,al,ar,m
  const double A       = al*B0(bl)+ar*B0(br);
  const double B0l     = B0(bl);
  const double B0r     = B0(br);
  const double B1l     = B1(bl);
  const double B1r     = B1(br);
  const double B2l     = B2(bl);
  const double B2r     = B2(br);
  const double dB0ldx  = dB0dx(bl);
  const double dB0rdx  = dB0dx(br);
  const double dB0ldx2 = dB0dx2(bl);
  const double dB0rdx2 = dB0dx2(br);

  // bl - bl
  I(0,0) = al*(dB0ldx2-al*dB0ldx*dB0ldx/A +
                            B2l/bl -2*B1l/(bl*bl)+2*B0l/pow(bl,3))/A;

  // bl - br
  I(0,1) = -al*ar*dB0ldx*dB0rdx/(A*A);
  I(1,0) = gsl_matrix_get(I,0,1);

  // bl - al
  I(0,2) = dB0ldx/A-al*B0l*dB0ldx/(A*A)-B1l/A;
  I(2,0) = gsl_matrix_get(I,0,2);

  // bl - ar
  I(0,3) = -al*B0r*dB0ldx/(A*A);
  I(3,0) = gsl_matrix_get(I,0,3);

  // br - br
  I(1,1) = ar*(dB0rdx2-ar*dB0rdx*dB0rdx/A +
                            B2r/br -2*B1r/(br*br)+2*B0r/pow(br,3))/A ;

  // br - al
  I(1,2) = -ar*B0l*dB0rdx/(A*A);
  I(2,1) = gsl_matrix_get(I,1,2);

  // br - ar
  I(1,3) = dB0rdx/A-ar*B0r*dB0rdx/(A*A)-B1r/A;
  I(3,1) = gsl_matrix_get(I,1,3);

  // al parameter
  // al - al
  I(2,2) = -B0l*B0l/(A*A)+(bl+1)*B0l/(al*A);

  // al - ar
  I(2,3) = -B0l*B0r/(A*A);
  I(3,2) = gsl_matrix_get(I,2,3);

  // ar parameter
  // ar - ar
  I(3,3) = -B0r*B0r/(A*A)+(br+1)*B0r/(ar*A);

  if(!O_munknown){

    // bl - m
    I(0,4) = (log(bl)-M_EULER)/(bl*A);
    I(4,0) = gsl_matrix_get(I,0,4);

    // br - m
    I(1,4) = -(log(br)-M_EULER)/(br*A);
    I(4,1) = gsl_matrix_get(I,1,4);

    // al - m
    I(2,4) = -bl/(al*A);
    I(4,2) = gsl_matrix_get(I,2,4);

    // ar - m
    I(3,4) = br/(ar*A);
    I(4,3) = gsl_matrix_get(I,3,4);

    const double dt1l = gsl_sf_gamma (2.-1/bl);
    const double dt1r = gsl_sf_gamma (2.-1/br);
    const double dt2l = pow(bl,1.-1./bl);
    const double dt2r = pow(br,1.-1./br);

    // m - m
    I(4,4) = (dt1l*dt2l/al+dt1r*dt2r/ar)/A;

    par = Rcpp::CharacterVector::create("bl", "br", "al", "ar", "m");

  }else{
    par = Rcpp::CharacterVector::create("bl", "br", "al", "ar");
  }

  // invert I; in J store temporary LU decomp.
  gsl_matrix_memcpy (J,I);
  gsl_linalg_LU_decomp (J,P,&signum);
  gsl_linalg_LU_invert (J,P,invI);


  // create vectors of standard errors
  Rcpp::NumericVector std_err;
  std_err = gsl_matrix_diagonal(invI);
  std_err = std_err/size;

  // add name of parameters
  Rcpp::DataFrame dt =
    Rcpp::DataFrame::create( Rcpp::Named("par")     = par
                            ,Rcpp::Named("std_err") = std_err
                            );

  // creates information matrix
  Rcpp::NumericMatrix infmatrix =
    Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(I));

  // creates inverse information matrix
  Rcpp::NumericMatrix inv_infmatrix =
    Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(invI));

  // add names for the matrices
  colnames(infmatrix)     = par;
  colnames(inv_infmatrix) = par;

  Rcpp::List ans =
    Rcpp::List::create(
                       Rcpp::Named("std_error")      = dt
                       ,Rcpp::Named("infmatrix")     = infmatrix
                       ,Rcpp::Named("inv_infmatrix") = inv_infmatrix
                       );

  // free allocated vector
  gsl_permutation_free(P);

  return ans;
}
