
#include "common.h"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

#include <gsl/gsl_linalg.h>

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>


//' Return the Fisher Information Matrix for the Subbotin Distribution
//'
//' Calculate the standard errors, the correlation, the
//' Fisher Information matrix and its inverse for a power exponential density
//' with given parameters
//'
//' @param size numeric - number of observations (Default: 01)
//' @param b numeric - the exponent b (Default: 2.0)
//' @param m numeric - the location parameter (Default: 0.0)
//' @param a numeric - the scale parameter (Default: 1.0)
//' @param O_munknown numeric - if true assumes m known
//' @return a list containing four elements:
//' * std_error - the standard error for the parameters
//' * cor_ab    - the correlation between parameters a and b
//' * infmatrix - the Fisher Information Matrix
//' * inv_infmatrix - the Inverse Fisher Information Matrix
//' @export
//' @md
// [[Rcpp::export]]
Rcpp::List subbofish(
                       unsigned size = 1
                      ,double b = 2.0
                      ,double m = 0.0
                      ,double a = 1.0
                      ,unsigned O_munknown = 0
                      ){

  // used only in the old code part
  //double g,psi1,psi2,var_b,var_a,var_m,
  //gsl_sf_result result;

  // correlation between a and b parameters
  double cor_ab;

  // dimension of the matrix
  size_t dim=(O_munknown?2:3);

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

  // the order is b,a,m
  const double dtmp1= log(b)+gsl_sf_psi(1+1./b);

  // b - b
  I(0,0) = (dtmp1*dtmp1+gsl_sf_psi_1(1+1/b)*(1+1/b)-1)/pow(b,3);

  // b - a
  I(0,1) = -dtmp1/(a*b);
  I(1,0) = gsl_matrix_get(I,0,1);

  // a - a
  I(1,1) = b/(a*a);

  if(!O_munknown){
  // b - m
    I(0,2) = 0;
    I(2,0) = gsl_matrix_get(I,0,2);

  // a - m
    I(1,2) = 0;
    I(2,1) = gsl_matrix_get(I,1,2);

  // m - m
    I(2,2) = (pow(b,-2/b+1)*gsl_sf_gamma(2-1/b))/
             (gsl_sf_gamma(1+1/b)*a*a);

    par = Rcpp::CharacterVector::create("b", "a", "m");

  }else{
    par = Rcpp::CharacterVector::create("b", "a");
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

  // create correlation between parameters a and b
  cor_ab = invI(0,1)/sqrt(invI(0,0)*invI(1,1));

  // OLD CODE: use the analytical inversion
  //       g=1.+1./b;
  //       gsl_sf_psi_1_e(g,&result);
  //       psi1=g*result.val-1.;
  //       psi2=log(b)+gsl_sf_psi(g);
  //       var_b = b*b*b/psi1;
  //       var_m = a*a*pow(b,2/b-1.)*gsl_sf_gamma(1.+1./b)/gsl_sf_gamma(2.-1./b);
  //       var_a = a*a*(1.+psi2*psi2/psi1)/b;
  //       cor_ab = (b*a*psi2/psi1);
  //       if(O_verbose)
  //    fprintf(stderr,"# sigma_b      sigma_a      sigma_m      corr_a_b\n");
  //       printf("% 12.4e % 12.4e % 12.4e % 12.4e\n",
  //         sqrt(var_b/size),sqrt(var_a/size),sqrt(var_m/size),cor_ab/sqrt(var_b*var_a));

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
                       ,Rcpp::Named("cor_ab")        = cor_ab
                       ,Rcpp::Named("infmatrix")     = infmatrix
                       ,Rcpp::Named("inv_infmatrix") = inv_infmatrix
                       );


  // free allocated vector
  gsl_permutation_free(P);

  return ans;
}
