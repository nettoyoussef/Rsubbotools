/*
  Probability Functions for several distributions (CDFs)
  likelihood maximization

  Copyright (C) 2007-2014 Giulio Bottazzi
  Copyright (C) 2020-2021 Elias Youssef Haddad Netto

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  (version 2) as published by the Free Software Foundation;

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/

#include "common.h"
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_integration.h>


// Auxiliar function
#include <algorithm>
#include <iostream>
#include <vector>

// //' Returns CDF from the Gamma Distribution
// //'
// //' The \code{pgamma_c} returns the Cumulative Distribution Function at point x
// //' for the Gamma distribution with parameters \eqn{a}, \eqn{b}.
// //'
// //' The gamma distribution is given by the function:
// //' \deqn{f(x) = \frac{1}{\Gamma(b)a^b}x^{b-1}e^{-x/a}, x > 0}
// //' where \eqn{b} is a shape parameter, \eqn{a} is a scale parameter, and
// //' \eqn{\Gamma} represents the gamma function.
// //' @param x (numeric) - value in the range \eqn{(0, \infty)} to evaluate the density.
// //' @param a (numeric) - scale parameter. Must be in the range \eqn{(0, \infty)}.
// //' @param b (numeric) - shape parameter. Must be in the range \eqn{(0, \infty)}.
// //' @return a vector containing the values for the densities.
// //' @export
// //' @md
// // [[Rcpp::export]]
// Rcpp::NumericVector pgamma_c(
//                              Rcpp::NumericVector x
//                              ,double a = 1
//                              ,double b = 2
//                              ){
//
//   // parameters check
//   if (a < 0 || b < 0){
//     Rcpp::stop("Parameters a and b must be greater than zero.");
//   }
//
//   unsigned int n, i;
//
//   // size of x vector
//   n = x.size();
//   Rcpp::NumericVector cdf(n);
//
//   for(i = 0; i < n; i++){
//     // x values' check
//     if(x[i] < 0){
//       Rprintf("x[%d] = %f", i, x[i]);
//       Rcpp::stop("All x values must be greater than zero.");
//     }
//     // calculates cdf
//     cdf[i] = gsl_sf_gamma_inc_P(b, x[i]/a);
//   }
//   return cdf;
// }



//' Returns CDF from EP Distribution
//'
//' The \code{ppower} returns the Cumulative Distribution Function at point x for
//' the Exponential Power distribution with parameters \eqn{a}, \eqn{b} and \eqn{m}.
//'
//' The Exponential Power distribution (EP) is given by the function:
//' \deqn{ f(a,b) = \frac{1}{2a\Gamma(1+1/b)}e^{-|(x-m)/a|^b}, -\infty < x < \infty }.
//' where \eqn{b} is a shape parameter, \eqn{a} is a scale parameter, \eqn{m}
//' is a location parameter and \eqn{\Gamma} represents the gamma function.
//' @param x (numeric) - value in the range \eqn{(-\infty, \infty)} to evaluate
//' the density.
//' @param m (numeric) - location parameter. Must be in the range
//' @param a (numeric) - scale parameter. Must be in the range \eqn{(0, \infty)}.
//' @param b (numeric) - shape parameter. Must be in the range \eqn{(0, \infty)}.
//' \eqn{(-\infty, \infty)}.
//' @return a vector containing the values for the densities.
//' @export
//' @md
// [[Rcpp::export]]
Rcpp::NumericVector ppower(
  Rcpp::NumericVector x
  ,double m = 0
  ,double a = 1
  ,double b = 2
){

  if (a < 0 || b < 0){
    Rcpp::stop("Parameters a and b must be greater than zero.");
  }
  unsigned int n, i;

  // size of x vector
  n = x.size();
  Rcpp::NumericVector cdf(n);
  double tmp;

  // calculates cdf
  for(i = 0; i < n; i++){
    tmp = (x[i] - m)/a;

    // debug
    //Rprintf("sgn: %d\n", sgn(x[i]-m));
    //Rprintf("gamma: %f\n", gsl_sf_gamma_inc_P(1./b, pow( fabs(tmp), b))/2);
    cdf[i] = 0.5 + sgn(x[i]-m)*gsl_sf_gamma_inc_P(1./b, pow( fabs(tmp), b))/2;
  }
  return cdf;
}


// density function for the SEP
double fSEP(double x, void *params){

  const double *X = (double *) params;

  const double m=X[0]; // location
  const double a=X[1]; // scale
  const double b=X[2]; // shape
  const double lambda=X[3]; // skewness

  const double z = (x-m)/a;
  const double w = sgn(z)*pow(fabs(z),b/2)*lambda*sqrt(2/b);
  const double c = 2*a*pow(b,1/b-1)*gsl_sf_gamma(1./b);

  return 2*gsl_cdf_ugaussian_P(w)*exp(-pow(fabs(z),b)/b)/c;
}


//' Returns CDF from the Skewed Exponential Power distribution
//'
//' The \code{psep} returns the Cumulative Distribution Function at point x for
//' the Skewed Exponential Power distribution with parameters \eqn{a}, \eqn{b}.
//'
//' The  SEP is a exponential power distribution controlled
//' by four parameters, with formula:
//' \deqn{ f(x; m, b, a, \lambda) = 2 \Phi(w) e^{-|z|^b/b}/(c)}
//' where:
//' \deqn{z = (x-m)/a}
//' \deqn{w = sign(z) |z|^{(b/2)} \lambda \sqrt{2/b}}
//' \deqn{c = 2 ab^{(1/b)-1} \Gamma(1/b)}
//' with \eqn{\Phi} the cumulative normal distribution with mean zero and variance
//' one. The CDF is calculated through numerical integration using the GSL suite.
//' @param x - vector with values to evaluate CDF.
//' @param m - the location parameter.
//' @param a - the scale parameter.
//' @param b - the shape parameter
//' @param lambda - the skewness parameter.
//' @export
//' @md
// [[Rcpp::export]]
Rcpp::NumericVector psep(
  Rcpp::NumericVector x
  ,double m      = 0
  ,double a      = 2
  ,double b      = 1
  ,double lambda = 0
){

  // need to declare a vector of parameters to use with GSL integration suite
  double par[4];
  par[0]=m;
  par[1]=a;
  par[2]=b;
  par[3]=lambda;

  size_t i,n;
  double error;

  // size of x vector
  n = x.size();
  Rcpp::NumericVector cdf(n);
  std::vector<size_t> indices(n);

  // saves data original position
  indices = rank_vector( Rcpp::as<std::vector<double>>(x) );

  // value check for R
  for(i=0; i <n; i++){
    if( x[i] == R_PosInf ){
      Rcpp::stop("Infinity values found on x[%d]", i);
    }
    if( x[i] == R_NegInf ){
      Rcpp::stop("-Infinity values found on x[%d]", i);
    }
  }

  /* allocate at most n subintervals */
  gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (10000);

  /* define the function to be integrated */
  gsl_function F;
  F.function = &fSEP;
  F.params = par;

  for(i=0;i<n;i++){
    gsl_integration_qagil(&F,x[i],1e-7,1e-7,10000,w,&cdf[i],&error);
    // the quadrature procedure can fail for extreme values of x, much above
    // the 99 percentile of the distribution
    if((x[i] > 3*m) && cdf[i] < 0.5){
      cdf[i] = 1;
    }
  }
  // free the workspace memory
  gsl_integration_workspace_free(w);

  return cdf;
}


//' Returns CDF from the Laplace Distribution
//'
//' The \code{plaplace} returns the Cumulative Distribution Function at point x
//' for the Laplace distribution with parameters \eqn{a} and \eqn{m}.
//'
//' The Laplace distribution is a distribution controlled
//' by two parameters, with formula:
//' \deqn{f(x;a,m) = \frac{1}{2a} e^{- \left| \frac{x-m}{a} \right| }}
//' where \eqn{a} is a scale parameter, and \eqn{m} is a location parameter.
//' @param x (numeric) - value in the range \eqn{(-\infty, \infty)} to evaluate
//' the density.
//' @param a (numeric) - scale parameter. Must be in the range \eqn{(0, \infty)}.
//' @param m (numeric) - location parameter.
//' @return a vector containing the values for the densities.
//' @export
//' @md
// [[Rcpp::export]]
Rcpp::NumericVector plaplace(
  Rcpp::NumericVector x
  ,double m = 0
  ,double a = 1
){

  if (a < 0){
    Rcpp::stop("Parameter a must be greater than zero.");
  }
  unsigned int n, i;

  // size of x vector
  n = x.size();
  Rcpp::NumericVector cdf(n);

  // calculates cdf
  for(i = 0; i < n; i++){
    if(x[i]>m){
      cdf[i] = 1-exp( (m-x[i])/a )/2;
    }
    else{
      cdf[i] = exp( (x[i]-m)/a )/2;
    }
  }
  return cdf;
}

//' Returns CDF from Asymmetric Laplace Distribution
//'
//' The \code{palaplace} returns the Cumulative Distribution Function at point x
//' for the Asymmetric Laplace distribution with parameters \eqn{a*} and \eqn{m}.
//'
//' The Asymmetric Laplace distribution is a distribution controlled
//' by three parameters, with formula:
//' \deqn{f(x;a_l,a_r,m) =
//' \frac{1}{A} e^{-|\frac{x-m}{a_l}| }, x < m
//' }
//' \deqn{f(x;a_l,a_r,m) =
//' \frac{1}{A} e^{-|\frac{x-m}{a_r}| }, x > m
//' }
//' with:
//' \deqn{A = a_l + a_r}
//' where \eqn{a*} are scale parameters, and \eqn{m} is a location parameter.
//' It is basically derived from the Asymmetric Exponential Power distribution
//' by setting \eqn{b_l = b_r = b}.
//' @param x (numeric) - value in the range \eqn{(-\infty, \infty)} to evaluate
//' the density.
//' @param m (numeric) - location parameter.
//' @param al,ar (numeric) - scale parameters. Must be in the range
//' \eqn{(0, \infty)}.
//' @return a vector containing the values for the densities.
//' @export
//' @md
// [[Rcpp::export]]
Rcpp::NumericVector palaplace(
  Rcpp::NumericVector x
  ,double m  = 0
  ,double al = 1
  ,double ar = 1
){

  if (al < 0 || ar < 0){
    Rcpp::stop("Parameters a must be greater than zero.");
  }
  unsigned int n, i;

  // size of x vector
  n = x.size();
  Rcpp::NumericVector cdf(n);

  // normalization constant
  double A = (al+ar);
  double tmp;

  // calculates cdf
  for(i = 0; i < n; i++){
    if(x[i]>m){
      tmp=(x[i]-m)/ar;
      cdf[i]=1-exp(-tmp)*ar/A;
    }
    else{
      tmp=(m-x[i])/al;
      cdf[i]=exp(-tmp)*al/A;
    }
  }
  return cdf;
}

//' Returns CDF from Subbotin Distribution
//'
//' The \code{psubbo} returns the Cumulative Distribution Function (CDF) from the
//' the Subbotin evaluated at \eqn{a} and return \eqn{z}, such that
//' \eqn{P(X < a) = z}.
//'
//' The Subbotin cumulative distribution function is given by:
//' \deqn{F(x;a,b,m) = 0.5 + 0.5 \text{sign}(x -m)P(x, 1/b)}
//' where \eqn{P} is the normalized incomplete gamma function:
//' \deqn{P(x, 1/b) = 1 - \frac{1}{\Gamma(1/b)} \int_{0}^{x} t^{1/b -1}e^{-t} }
//' and \eqn{a} is a scale parameter, \eqn{b} controls the tails (lower values
//' represent fatter tails), and \eqn{m} is a location parameter.
//' @param x (numeric) - value in the range \eqn{(-\infty, \infty)} to evaluate
//' the density.
//' @param m (numeric) - location parameter.
//' @param a (numeric) - scale parameter. Must be in the range \eqn{(0, \infty)}.
//' @param b (numeric) - shape parameter. Must be in the range \eqn{(0, \infty)}.
//' @return a vector containing the values for the densities.
//' @export
//' @md
// [[Rcpp::export]]
Rcpp::NumericVector psubbo(
  Rcpp::NumericVector x
  ,double m = 0
  ,double a = 1
  ,double b = 2
){

  if (a < 0 || b < 0){
    Rcpp::stop("Parameters a and b must be greater than zero.");
  }
  unsigned int n, i;
  double tmp;

  // size of x vector
  n = x.size();
  Rcpp::NumericVector cdf(n);

  // calculates density
  for(i = 0; i < n; i++){
    tmp = pow(fabs(x[i]-m)/a, b)/b;
    cdf[i]= .5+.5*(sgn(x[i]-m))*gsl_sf_gamma_inc_P(1./b, tmp);
  }
  return cdf;
}


//' Returns CDF from the AEP Distribution
//'
//' The \code{pasubbo} returns the Cumulative Distribution Function at point x
//' for the AEP distribution with parameters \eqn{a*}, \eqn{b*}, \eqn{m}.
//'
//' The AEP is a exponential power distribution controlled
//' by five parameters, with formula:
//' \deqn{ f(x;a_l,a_r,b_l,b_r,m) =
//' \frac{1}{A} e^{- \frac{1}{b_l} |\frac{x-m}{a_l}|^{b_l} }, x < m
//' }
//' \deqn{ f(x;a_l,a_r,b_l,b_r,m) =
//' \frac{1}{A} e^{- \frac{1}{b_r} |\frac{x-m}{a_r}|^{b_r} }, x > m
//' }
//' with:
//' \deqn{A = a_lb_l^{1/b_l}\Gamma(1+1/b_l) + a_rb_r^{1/b_r}\Gamma(1+1/b_r)}
//' where \eqn{l} and \eqn{r} represent left and right tails, \eqn{a*} are
//' scale parameters, \eqn{b*} control the tails (lower values represent
//' fatter tails), and \eqn{m} is a location parameter.
//'
//' @param x (numeric) - value in the range \eqn{(-\infty, \infty)} to evaluate
//' the density.
//' @param m (numeric) - location parameter.
//' @param al,ar (numeric) - scale parameters. Must be in the range
//' \eqn{(0, \infty)}.
//' @param bl,br (numeric) - shape parameters. Must be in the range
//' \eqn{(0, \infty)}.
//' @return a vector containing the values for the densities.
//' @export
//' @md
// [[Rcpp::export]]
Rcpp::NumericVector pasubbo(
  Rcpp::NumericVector x
  ,double m  = 0
  ,double al = 1
  ,double ar = 1
  ,double bl = 2
  ,double br = 2
){

  if (al < 0 || ar < 0 || bl < 0 || br < 0){
    Rcpp::stop("Parameters al, ar, bl, br must be greater than zero.");
  }

  unsigned int n, i;
  double tmp;

  // normalization constant
  const double Al=al*gsl_sf_gamma(1./bl+1.)*pow(bl,1./bl);
  const double Ar=ar*gsl_sf_gamma(1./br+1.)*pow(br,1./br);
  const double Asum=Al+Ar;

  // size of x vector
  n = x.size();
  Rcpp::NumericVector cdf(n);

  // calculates cdf
  for(i=0; i<n; i++){
    if(x[i]>m){
      tmp=pow((x[i]-m)/ar,br)/br;
      cdf[i]=(Al+Ar*gsl_sf_gamma_inc_P(1./br,tmp))/Asum;
     }
    else{
      tmp=pow((m-x[i])/al,bl)/bl;
      cdf[i]=Al*gsl_sf_gamma_inc_Q(1./bl,tmp)/Asum;
    }
  }
  return cdf;
}
