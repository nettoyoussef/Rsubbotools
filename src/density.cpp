#include "common.h"
#include <gsl/gsl_cdf.h>

// //' Returns density from Gamma Distribution
// //'
// //' The \code{dgamma_c} returns the density at point x for the
// //' Gamma distribution with parameters \eqn{a}, \eqn{b}.
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
// Rcpp::NumericVector dgamma_c(
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
//   Rcpp::NumericVector density(n);
//
//   // x values' check
//     for(i = 0; i < n; i++){
//       if(x[i] < 0){
//         Rprintf("x[%d] = %f", i, x[i]);
//         Rcpp::stop("All x values must be greater than zero.");
//       }
//
//   }
//
//   // normalization constant
//   double A = pow(a, b)*gsl_sf_gamma(b);
//
//   // calculates density
//   for(i = 0; i < n; i++){
//     density[i] = pow(x[i], b-1)*exp(-x[i]/a)/A;
//   }
//
//   return density;
// }


//' Returns density from EP Distribution
//'
//' The \code{dpower} returns the density at point x for the
//' Exponential Power distribution with parameters \eqn{a}, \eqn{b} and \eqn{m}.
//'
//' The Exponential Power distribution (EP) is given by the function:
//' \deqn{ f(a,b) = \frac{1}{2a\Gamma(1+1/b)}e^{-|(x-m)/a|^b}, -\infty < x < \infty }.
//' where \eqn{b} is a shape parameter, \eqn{a} is a scale parameter, \eqn{m}
//' is a location parameter and \eqn{\Gamma} represents the gamma function.
//' @param x (numeric) - value in the range \eqn{(-\infty, \infty)} to evaluate
//' the density.
//' @param m (numeric) - location parameter. Must be in the range
//' \eqn{(-\infty, \infty)}.
//' @param a (numeric) - scale parameter. Must be in the range \eqn{(0, \infty)}.
//' @param b (numeric) - shape parameter. Must be in the range \eqn{(0, \infty)}.
//' @return a vector containing the values for the densities.
//' @export
//' @md
// [[Rcpp::export]]
Rcpp::NumericVector dpower(
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
  Rcpp::NumericVector density(n);

  // normalization constant
  double A = 2*a*gsl_sf_gamma(1.+1./b);

  // calculates density
  for(i = 0; i < n; i++){
    density[i] = exp(-pow(fabs(x[i]-m)/a, b))/A;
  }
  return density;
}

//' Returns density from Skewed Exponential Power distribution
//'
//' The \code{dsep} returns the density at point x for the
//' Gamma distribution with parameters \eqn{a}, \eqn{b}.
//'
//' The  SEP is a exponential power distribution controlled
//' by four parameters, with formula:
//' \deqn{ f(x; \mu, \alpha, \lambda, \sigma) =
//' 2 \Phi(w) e^{-|z|^\alpha/\alpha}/ ( \sigma C)}
//' where:
//' \deqn{z = (x-\mu)/\sigma}
//' \deqn{w = sign(z) |z|^{(\alpha/2)} \lambda \sqrt{2/\alpha}}
//' \deqn{C = 2 \alpha^{(1/\alpha-1)} \Gamma(1/\alpha)}
//' with \eqn{\Phi} the cumulative normal distribution with mean zero and variance
//' one.
//' @param x (numeric) - value in the range \eqn{(-\infty, \infty)} to evaluate
//' the density.
//' @param m (numeric) - location parameter. Must be in the range
//' \eqn{(-\infty, \infty)}.
//' @param a (numeric) - scale parameter. Must be in the range \eqn{(0, \infty)}.
//' @param b (numeric) - shape parameter. Must be in the range \eqn{(0, \infty)}.
//' @param lambda (numeric) - skewness parameter. Must be in the range \eqn{(-\infty, \infty)}.
//' @return a vector containing the values for the densities.
//' @export
//' @md
// [[Rcpp::export]]
Rcpp::NumericVector dsep(
  Rcpp::NumericVector x
  ,double m      = 0
  ,double a      = 1
  ,double b      = 1
  ,double lambda = 1
){

  if (a < 0 || b < 0){
    Rcpp::stop("Parameters a and b must be greater than zero.");
  }

  unsigned int n, i;

  // size of x vector
  n = x.size();
  Rcpp::NumericVector density(n);
  double z;
  double w;
  double c;

  // calculates density
  for(i = 0; i < n; i++){
    z = (x[i]-m)/a;
    w = sgn(z)*pow(fabs(z),b/2)*lambda*sqrt(2/b);
    c = 2*a*pow(b,1/b-1)*gsl_sf_gamma(1./b);
    density[i] = 2*gsl_cdf_ugaussian_P(w)*exp(-pow(fabs(z),b)/b)/c;
    //dpower(x[i], b*pow(a, 1./a -1), a);
  }

  return density;
}


//' Returns density from Laplace Distribution
//'
//' The \code{dlaplace} returns the density at point x for the
//' Laplace distribution with parameters \eqn{a} and \eqn{m}.
//'
//' The Laplace distribution is a distribution controlled
//' by two parameters, with formula:
//' \deqn{f(x;a,m) = \frac{1}{2a} e^{- \left| \frac{x-m}{a} \right| }}
//' where \eqn{a} is a scale parameter, and \eqn{m} is a location parameter.
//' @param x (numeric) - value in the range \eqn{(-\infty, \infty)} to evaluate
//' the density.
//' @param m (numeric) - location parameter.
//' @param a (numeric) - scale parameter. Must be in the range \eqn{(0, \infty)}.
//' @return a vector containing the values for the densities.
//' @export
//' @md
// [[Rcpp::export]]
Rcpp::NumericVector dlaplace(
  Rcpp::NumericVector x
  ,double m = 0
  ,double a = 1
){

  if (a <= 0){
    Rcpp::stop("Parameter a must be greater than zero.");
  }
  unsigned int n, i;

  // size of x vector
  n = x.size();
  Rcpp::NumericVector density(n);

  // calculates density
  for(i = 0; i < n; i++){
    density[i] = exp(-fabs(x[i]-m)/a)/(2*a);
  }
  return density;

}


//' Returns density from Asymmetric Laplace Distribution
//'
//' The \code{dalaplace} returns the density at point x for the
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
//' @param x (numeric) - value in the range \eqn{(-\infty, \infty)} to evaluate
//' the density.
//' @param m (numeric) - location parameter.
//' @param al,ar (numeric) - scale parameters. Must be in the range
//' \eqn{(0, \infty)}.
//' @return a vector containing the values for the densities.
//' @export
//' @md
// [[Rcpp::export]]
Rcpp::NumericVector dalaplace(
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
  Rcpp::NumericVector density(n);

  // normalization constant
  double A = (al+ar);

  // calculates density
  for(i = 0; i < n; i++){
    if(x[i] < m){
      density[i] = exp(-fabs(x[i]-m)/al)/A;
    }else{
      density[i] = exp(-fabs(x[i]-m)/ar)/A;
    }
  }
  return density;
}

//' Returns density from Subbotin Distribution
//'
//' The \code{dsubbo} returns the density at point x for the
//' Subbotin distribution with parameters \eqn{a}, \eqn{b}, \eqn{m}.
//'
//' The Subbotin distribution is a exponential power distribution controlled
//' by three parameters, with formula:
//' \deqn{f(x;a,b,m) = \frac{1}{A} e^{-\frac{1}{b} |\frac{x-m}{a}|^b}}
//' with:
//' \deqn{A = 2ab^{1/b}\Gamma(1+1/b)}
//' where \eqn{a} is a scale parameter, \eqn{b} controls the tails (lower values
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
Rcpp::NumericVector dsubbo(
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
  Rcpp::NumericVector density(n);

  // normalization constant
  double A = (2*a*pow(b, 1./b)*gsl_sf_gamma(1.+1./b));

  // calculates density
  for(i = 0; i < n; i++){
    density[i] = exp(-pow(fabs(x[i]-m)/a, b)/b)/A;
  }

  return density;
}

//' Returns density from the AEP Distribution
//'
//' The \code{dasubbo} returns the density at point x for the
//' AEP distribution with parameters \eqn{a*}, \eqn{b*}, \eqn{m}. Notice
//' that this function can generate RNGs for both the \code{subboafit} and
//' \code{subbolafit} routines.
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
Rcpp::NumericVector dasubbo(
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

  // normalization constant
  double A = al*pow(bl,1/bl)*gsl_sf_gamma(1./bl+1.)
           + ar*pow(br,1/br)*gsl_sf_gamma(1./br+1.);

  // size of x vector
  n = x.size();
  Rcpp::NumericVector density(n);

  // calculates density
  for(i = 0; i < n; i++){
    if(x[i] < m){
      density[i] = exp(-pow(fabs(x[i]-m)/al, bl)/bl)/A;
    }else{
      density[i] = exp(-pow(fabs(x[i]-m)/ar, br)/br)/A;
    }
  }
  return density;
}
