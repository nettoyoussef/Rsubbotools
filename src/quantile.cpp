
// Cumulative Distribution Functions Functions

#include "common.h"
#include "cdf.h"
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_integration.h>

// Auxiliar function
#include <algorithm>
#include <iostream>
#include <vector>


//' Returns quantile from EP Distribution
//'
//' The \code{qpower} returns the density at point x for the
//' Exponential Power distribution with parameters \eqn{a}, \eqn{b} and \eqn{m}.
//'
//' The Exponential Power distribution (EP) is given by the function:
//' \deqn{ f(a,b) = \frac{1}{2a\Gamma(1+1/b)}e^{-|(x-m)/a|^b}, -\infty < x < \infty }.
//' where \eqn{b} is a shape parameter, \eqn{a} is a scale parameter, \eqn{m}
//' is a location parameter and \eqn{\Gamma} represents the gamma function.
//' Copyright (C) 2021 Elias Haddad
//'
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
Rcpp::NumericVector qpower(
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
  Rcpp::NumericVector quantile(n);

  // calculates cdf
  for(i = 0; i < n; i++){
    if(x[i] < 0 || x[i] > 1){
      Rprintf("x[%d]= %f", i, x[i]);
      Rcpp::stop("x range must be in [0, 1].");
    }
    quantile[i] = (
      m + 
      a 
      * sgn(x[i]-0.5)
      * pow(
         inv_inc_lower_gamma(1./b, 2*gsl_sf_gamma(1./b)*fabs(x[i] - 0.5))
         , 1./b
      )
     );
  }
  return quantile;
}

double diff_psep(newton_args x){

  double m, b, a, lambda, result;
  Rcpp::NumericVector y(1);
  y[0]   = x.x_guess;
  m      = x.func_args[0];
  a      = x.func_args[1];
  b      = x.func_args[2];
  lambda = x.func_args[3];

  // since psep is vectorized, and this function isn't
  // I have to always get the first element
  Rcpp::NumericVector tmp;
  tmp = psep(y, m, a, b, lambda);
  //Rprintf("tmp: %.4f\n", tmp[0]);
  //Rprintf("x_guess: %.4f\n", x.x_guess);
  //Rprintf("m: %.4f\n", m);
  //Rprintf("b: %.4f\n", b);
  //Rprintf("a: %.4f\n", a);
  //Rprintf("lambda: %.4f\n", lambda);
  result = tmp[0] - x.x;

  return result;
}

double pdf_sep(newton_args x){

  double m=x.func_args[0]; // location
  double a=x.func_args[1]; // scale
  double b=x.func_args[2]; // shape
  double lambda=x.func_args[3]; // skewness

  double z = (x.x_guess-m)/a;
  double w = sgn(z)*pow(fabs(z),b/2)*lambda*sqrt(2/b);
  double c = 2*a*pow(b,1/b-1)*gsl_sf_gamma(1./b);

  return 2*gsl_cdf_ugaussian_P(w)*exp(-pow(fabs(z),b)/b)/c;

}

//' Returns quantile from the Skewed Exponential Power distribution
//'
//' The \code{qsep} returns the Cumulative Distribution Function at point x for
//' the Skewed Exponential Power distribution with parameters \eqn{a}, \eqn{b}.
//'
//' The  SEP is a exponential power distribution controlled
//' by four parameters, with formula:
//' \deqn{ f(x; \mu, a, \lambda, b) =
//' 2 \Phi(w) e^{-|z|^a/a}/ ( b C)}
//' where:
//' \deqn{z = (x-\mu)/b}
//' \deqn{w = sign(z) |z|^{(a/2)} \lambda \sqrt{2/a}}
//' \deqn{C = 2 a^{(1/a-1)} \Gamma(1/a)}
//' with \eqn{\Phi} the cumulative normal distribution with mean zero and variance
//' one. The CDF is calculated through numerical integration using the GSL suite.
//' Copyright (C) 2007-2014 Giulio Bottazzi
//' Copyright (C) 2020-2021 Elias Haddad
//' @param x (numeric) - vector with values to evaluate CDF.
//' @param m (numeric) - the location parameter.
//' @param a (numeric) - the scale parameter.
//' @param b (numeric) - the shape parameter
//' @param lambda (numeric) - the skewness parameter.
//' @param method (numeric) - If 0, uses the Newton-Raphson procedure for optimization. If 1, uses Steffensen.
//' @param step_size (numeric) - the size of the step in the numerical optimization (gradient descent). Default is 1e-4.
//' @param tol (numeric) - error tolerance (default is 1e-10).
//' @param max_iter (numeric) - maximum number of iterations for the optimization procedure (default is 100).
//' @param verb (numeric) - verbosity level of the process (default 0).
//' @return a vector containing the values for the densities.
//' @export
//' @md
// [[Rcpp::export]]
Rcpp::NumericVector qsep(
  Rcpp::NumericVector x
  ,double m         = 0
  ,double a         = 2
  ,double b         = 1
  ,double lambda    = 0
  ,int method       = 0
  ,double step_size = 1e-4
  ,double tol       = 1e-10
  ,int max_iter     = 100
  ,int verb         = 0
){

  int i,n;
  n = x.size();
  Rcpp::NumericVector quantile(n);
  Rcpp::NumericVector func_args(4);
  func_args[0] = m;
  func_args[1] = a;
  func_args[2] = b;
  func_args[3] = lambda;

  newton_args z;
  z.func_args = func_args;
  z.f         = diff_psep;
  z.step_size = step_size;
  z.tol       = tol;
  z.max_iter  = max_iter;
  z.verb      = verb;
  z.dfdx      = pdf_sep;

  // Newton's method
  // Default, since it was faster on tests:
  // Unit: microseconds
  //                expr      min        lq      mean    median       uq      max  neval
  // qsep(b, method = 0)  584.814  620.9425  666.1068  640.7215  681.827 3454.803  10000
  // qsep(b, method = 1) 1000.355 1063.8695 1139.8203 1105.0830 1162.220 3383.893  10000
  if(method == 0){
    // value check for R
    for(i=0; i <n; i++){
      if(x[i] == 0){
      quantile[i] = R_NegInf;
      }else if(x[i] == 1){
      quantile[i] = R_PosInf;
      }else{
      // initial guess - starts by the average
      z.x_guess   = m;
      z.x         = x[i];
      quantile[i] = newton_c(z);
      }
    }
  }

  // Steffensen's method
  if(method == 1){
    // value check for R
    for(i=0; i <n; i++){
      if(x[i] == 0){
      quantile[i] = R_NegInf;
      }else if(x[i] == 1){
      quantile[i] = R_PosInf;
      }else{
      // initial guess - starts by the average
      z.x_guess   = m;
      z.x         = x[i];
      quantile[i] = steffensen_c(z);
      }
    }
  }

  return quantile;
}

//' Returns quantile from Laplace Distribution
//'
//' The \code{qlaplace} returns the density at point x for the
//' Laplace distribution with parameters \eqn{a} and \eqn{m}.
//'
//' The Laplace distribution is a distribution controlled
//' by two parameters, with formula:
//' \deqn{f(x;a,m) = \frac{1}{2a} e^{- \left| \frac{x-m}{a} \right| }}
//' where \eqn{a} is a scale parameter, and \eqn{m} is a location parameter.
//' Copyright (C) 2021 Elias Haddad
//'
//' @param x (numeric) - value in the range \eqn{(-\infty, \infty)} to evaluate
//' the density.
//' @param m (numeric) - location parameter.
//' @param a (numeric) - scale parameter. Must be in the range \eqn{(0, \infty)}.
//' @return a vector containing the values for the densities.
//' @export
//' @md
// [[Rcpp::export]]
Rcpp::NumericVector qlaplace(
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
  Rcpp::NumericVector quantile(n);

  // calculates cdf
  for(i = 0; i < n; i++){
    quantile[i] = m - a*sgn(x[i]-0.5)*log(1 - 2*fabs(x[i]-0.5));
  }
  return quantile;
}

//' Returns quantile from Asymmetric Laplace Distribution
//'
//' The \code{qalaplace} returns the density at point x for the
//' Asymmetric Laplace distribution with parameters \eqn{a*} and \eqn{m}.
//'
//' The Asymmetric Laplace distribution is a distribution controlled
//' by three parameters, with formula:
//' \deqn{f(x;a_l,a_r,m) =
//' \begin{cases}
//' \frac{1}{A} e^{-|\frac{x-m}{a_l}| }, & x < m
//' \frac{1}{A} e^{-|\frac{x-m}{a_r}| }, & x > m
//' \end{cases}}
//' with:
//' \deqn{A = a_l + a_r}
//' where \eqn{a*} are scale parameters, and \eqn{m} is a location parameter.
//' It is basically derived from the Asymmetric Exponential Power distribution
//' by setting \eqn{b_l = b_r = b}.
//' Copyright (C) 2021 Elias Haddad
//'
//' @param x (numeric) - value in the range \eqn{(-\infty, \infty)} to evaluate
//' the density.
//' @param m (numeric) - location parameter.
//' @param al,ar (numeric) - scale parameters. Must be in the range
//' \eqn{(0, \infty)}.
//' @return a vector containing the values for the densities.
//' @export
//' @md
// [[Rcpp::export]]
Rcpp::NumericVector qalaplace(
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
  Rcpp::NumericVector quantile(n);

  double A = al+ar;

  // calculates quantile
  for(i = 0; i < n; i++){
    if(x[i] < al/A){
      quantile[i] = m + al*log(A*x[i]/al);
    }
    else{
      quantile[i] = m - ar*log(A*(1-x[i])/ar);
    }
  }
  return quantile;
}


//' Returns CDF from Subbotin Distribution
//'
//' The \code{qsubbo} returns the Cumulative Distribution Function (CDF) from the
//' the Subbotin evaluated at \eqn{a} and return \eqn{z}, such that
//' \eqn{P(X < a) = z}.
//'
//' The Subbotin cumulative distribution function is given by:
//' \deqn{F(x;a,b,m) = 0.5 + 0.5 \text{sign}(x -m)P(x, 1/b)}
//' where \eqn{P} is the normalized incomplete gamma function:
//' \deqn{P(x, 1/b) = 1 - \frac{1}{\Gamma(1/b)} \int_{0}^{x} t^{1/b -1}e^{-t} }
//' and \eqn{a} is a scale parameter, \eqn{b} controls the tails (lower values
//' represent fatter tails), and \eqn{m} is a location parameter.
//' Copyright (C) 2021 Elias Haddad
//'
//' @param x (numeric) - value in the range \eqn{(-\infty, \infty)} to evaluate
//' the density.
//' @param m (numeric) - location parameter.
//' @param a (numeric) - scale parameter. Must be in the range \eqn{(0, \infty)}.
//' @param b (numeric) - shape parameter. Must be in the range \eqn{(0, \infty)}.
//' @return a vector containing the values for the densities.
//' @export
//' @md
// [[Rcpp::export]]
Rcpp::NumericVector qsubbo(
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
  Rcpp::NumericVector quantile(n);

  // calculates density
  for(i = 0; i < n; i++){
    tmp = 2*fabs(x[i]-0.5)*gsl_sf_gamma(1./b);
    quantile[i]= m + a*sgn(x[i]-0.5)*pow(b*inv_inc_lower_gamma(1./b, tmp), 1./b);
  }
  return quantile;
}


//' Returns CDF from the AEP Distribution
//'
//' The \code{qasubbo} returns the density at point x for the
//' AEP distribution with parameters \eqn{a*}, \eqn{b*}, \eqn{m}.
//'
//' The AEP is a exponential power distribution controlled
//' by five parameters, with formula:
//' \deqn{ f(x;a_l,a_r,b_l,b_r,m) =
//' \begin{cases}
//' \frac{1}{A} e^{- \frac{1}{b_l} |\frac{x-m}{a_l}|^{b_l} }, & x < m
//' \frac{1}{A} e^{- \frac{1}{b_r} |\frac{x-m}{a_r}|^{b_r} }, & x > m
//' \end{cases} }
//' with:
//' \deqn{A = a_lb_l^{1/b_l}\Gamma(1+1/b_l) + a_rb_r^{1/b_r}\Gamma(1+1/b_r)}
//' where \eqn{l} and \eqn{r} represent left and right tails, \eqn{a*} are
//' scale parameters, \eqn{b*} control the tails (lower values represent
//' fatter tails), and \eqn{m} is a location parameter.
//' Copyright (C) 2020-2021 Elias Haddad
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
Rcpp::NumericVector qasubbo(
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
  Rcpp::NumericVector quantile(n);

  // calculates cdf
  for(i=0; i<n; i++){
    if(x[i] < Al/Asum ){
      tmp= inv_inc_upper_gamma(1./bl, gsl_sf_gamma(1./bl)*Asum*x[i]/Al);
      quantile[i]= m - al*pow(bl*tmp, 1./bl);
     }
    else{
      tmp= inv_inc_lower_gamma(1./br, gsl_sf_gamma(1./br)*Asum*(x[i] - Al/Asum)/Ar);
      quantile[i]= m + ar*pow(br*tmp, 1./br);
    }
  }
  return quantile;
}
