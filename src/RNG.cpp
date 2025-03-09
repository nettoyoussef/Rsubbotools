#include "common.h"

// tabulated values for the height of the Ziggurat levels
static const double ytab[128] = {
  1              , 0.963598623011  , 0.936280813353  , 0.913041104253 ,
  0.892278506696 , 0.873239356919  , 0.855496407634  , 0.838778928349 ,
  0.822902083699 , 0.807732738234  , 0.793171045519  , 0.779139726505 ,
  0.765577436082 , 0.752434456248  , 0.739669787677  , 0.727249120285 ,
  0.715143377413 , 0.703327646455  , 0.691780377035  , 0.68048276891  ,
  0.669418297233 , 0.65857233912   , 0.647931876189  , 0.637485254896 ,
  0.62722199145  , 0.617132611532  , 0.607208517467  , 0.597441877296 ,
  0.587825531465 , 0.578352913803  , 0.569017984198  , 0.559815170911 ,
  0.550739320877 , 0.541785656682  , 0.532949739145  , 0.524227434628 ,
  0.515614886373 , 0.507108489253  , 0.498704867478  , 0.490400854812 ,
  0.482193476986 , 0.47407993601   , 0.466057596125  , 0.458123971214 ,
  0.450276713467 , 0.442513603171  , 0.434832539473  , 0.427231532022 ,
  0.419708693379 , 0.41226223212   , 0.404890446548  , 0.397591718955 ,
  0.390364510382 , 0.383207355816  , 0.376118859788  , 0.369097692334 ,
  0.362142585282 , 0.355252328834  , 0.348425768415  , 0.341661801776 ,
  0.334959376311 , 0.328317486588  , 0.321735172063  , 0.31521151497  ,
  0.308745638367 , 0.302336704338  , 0.29598391232   , 0.289686497571 ,
  0.283443729739 , 0.27725491156   , 0.271119377649  , 0.265036493387 ,
  0.259005653912 , 0.253026283183  , 0.247097833139  , 0.241219782932 ,
  0.235391638239 , 0.229612930649  , 0.223883217122  , 0.218202079518 ,
  0.212569124201 , 0.206983981709  , 0.201446306496  , 0.195955776745 ,
  0.190512094256 , 0.185114984406  , 0.179764196185  , 0.174459502324 ,
  0.169200699492 , 0.1639876086    , 0.158820075195  , 0.153697969964 ,
  0.148621189348 , 0.143589656295  , 0.138603321143  , 0.133662162669 ,
  0.128766189309 , 0.123915440582  , 0.119109988745  , 0.114349940703 ,
  0.10963544023  , 0.104966670533  , 0.100343857232  , 0.0957672718266,
  0.0912372357329, 0.0867541250127 , 0.082318375932  , 0.0779304915295,
  0.0735910494266, 0.0693007111742 , 0.065060233529  , 0.0608704821745,
  0.056732448584 , 0.05264727098   , 0.0486162607163 , 0.0446409359769,
  0.0407230655415, 0.0368647267386 , 0.0330683839378 , 0.0293369977411,
  0.0256741818288, 0.0220844372634 , 0.0185735200577 , 0.0151490552854,
  0.0118216532614, 0.00860719483079, 0.00553245272614, 0.00265435214565
};

// tabulated values for 2^24 times x[i]/x[i+1],
// used to accept for U*x[i+1]<=x[i] without any floating point operations
static const unsigned long ktab[128] = {
  0       , 12590644, 14272653, 14988939,
  15384584, 15635009, 15807561, 15933577,
  16029594, 16105155, 16166147, 16216399,
  16258508, 16294295, 16325078, 16351831,
  16375291, 16396026, 16414479, 16431002,
  16445880, 16459343, 16471578, 16482744,
  16492970, 16502368, 16511031, 16519039,
  16526459, 16533352, 16539769, 16545755,
  16551348, 16556584, 16561493, 16566101,
  16570433, 16574511, 16578353, 16581977,
  16585398, 16588629, 16591685, 16594575,
  16597311, 16599901, 16602354, 16604679,
  16606881, 16608968, 16610945, 16612818,
  16614592, 16616272, 16617861, 16619363,
  16620782, 16622121, 16623383, 16624570,
  16625685, 16626730, 16627708, 16628619,
  16629465, 16630248, 16630969, 16631628,
  16632228, 16632768, 16633248, 16633671,
  16634034, 16634340, 16634586, 16634774,
  16634903, 16634972, 16634980, 16634926,
  16634810, 16634628, 16634381, 16634066,
  16633680, 16633222, 16632688, 16632075,
  16631380, 16630598, 16629726, 16628757,
  16627686, 16626507, 16625212, 16623794,
  16622243, 16620548, 16618698, 16616679,
  16614476, 16612071, 16609444, 16606571,
  16603425, 16599973, 16596178, 16591995,
  16587369, 16582237, 16576520, 16570120,
  16562917, 16554758, 16545450, 16534739,
  16522287, 16507638, 16490152, 16468907,
  16442518, 16408804, 16364095, 16301683,
  16207738, 16047994, 15704248, 15472926
};

// tabulated values of 2^{-24}*x[i]
static const double wtab[128] = {
  1.62318314817e-08, 2.16291505214e-08, 2.54246305087e-08, 2.84579525938e-08,
  3.10340022482e-08, 3.33011726243e-08, 3.53439060345e-08, 3.72152672658e-08,
  3.8950989572e-08 , 4.05763964764e-08, 4.21101548915e-08, 4.35664624904e-08,
  4.49563968336e-08, 4.62887864029e-08, 4.75707945735e-08, 4.88083237257e-08,
  5.00063025384e-08, 5.11688950428e-08, 5.22996558616e-08, 5.34016475624e-08,
  5.44775307871e-08, 5.55296344581e-08, 5.65600111659e-08, 5.75704813695e-08,
  5.85626690412e-08, 5.95380306862e-08, 6.04978791776e-08, 6.14434034901e-08,
  6.23756851626e-08, 6.32957121259e-08, 6.42043903937e-08, 6.51025540077e-08,
  6.59909735447e-08, 6.68703634341e-08, 6.77413882848e-08, 6.8604668381e-08 ,
  6.94607844804e-08, 7.03102820203e-08, 7.11536748229e-08, 7.1991448372e-08 ,
  7.2824062723e-08 , 7.36519550992e-08, 7.44755422158e-08, 7.52952223703e-08,
  7.61113773308e-08, 7.69243740467e-08, 7.77345662086e-08, 7.85422956743e-08,
  7.93478937793e-08, 8.01516825471e-08, 8.09539758128e-08, 8.17550802699e-08,
  8.25552964535e-08, 8.33549196661e-08, 8.41542408569e-08, 8.49535474601e-08,
  8.57531242006e-08, 8.65532538723e-08, 8.73542180955e-08, 8.8156298059e-08 ,
  8.89597752521e-08, 8.97649321908e-08, 9.05720531451e-08, 9.138142487e-08  ,
  9.21933373471e-08, 9.30080845407e-08, 9.38259651738e-08, 9.46472835298e-08,
  9.54723502847e-08, 9.63014833769e-08, 9.71350089201e-08, 9.79732621669e-08,
  9.88165885297e-08, 9.96653446693e-08, 1.00519899658e-07, 1.0138063623e-07 ,
  1.02247952126e-07, 1.03122261554e-07, 1.04003996769e-07, 1.04893609795e-07,
  1.05791574313e-07, 1.06698387725e-07, 1.07614573423e-07, 1.08540683296e-07,
  1.09477300508e-07, 1.1042504257e-07 , 1.11384564771e-07, 1.12356564007e-07,
  1.13341783071e-07, 1.14341015475e-07, 1.15355110887e-07, 1.16384981291e-07,
  1.17431607977e-07, 1.18496049514e-07, 1.19579450872e-07, 1.20683053909e-07,
  1.21808209468e-07, 1.2295639141e-07 , 1.24129212952e-07, 1.25328445797e-07,
  1.26556042658e-07, 1.27814163916e-07, 1.29105209375e-07, 1.30431856341e-07,
  1.31797105598e-07, 1.3320433736e-07 , 1.34657379914e-07, 1.36160594606e-07,
  1.37718982103e-07, 1.39338316679e-07, 1.41025317971e-07, 1.42787873535e-07,
  1.44635331499e-07, 1.4657889173e-07 , 1.48632138436e-07, 1.50811780719e-07,
  1.53138707402e-07, 1.55639532047e-07, 1.58348931426e-07, 1.61313325908e-07,
  1.64596952856e-07, 1.68292495203e-07, 1.72541128694e-07, 1.77574279496e-07,
  1.83813550477e-07, 1.92166040885e-07, 2.05295471952e-07, 2.22600839893e-07
};

// The Gamma distribution of order a>0 is defined by:
// p(x) dx = {1 / \Gamma(a) b^a } x^{a-1} e^{-x/b} dx
// for x>0.  If X and Y are independent gamma-distributed random
// variables of order a1 and a2 with the same scale parameter b, then
// X+Y has gamma distribution of order a1+a2.
// The algorithms below are from Knuth, vol 2, 2nd ed, p. 129.


//double gsl_ran_gaussian_ziggurat(const gsl_rng * r, double sigma)
//double gsl_ran_gaussian_ratio_method(const gsl_rng * r, double sigma)
//
// This function computes a Gaussian random variate using the alternative
// Marsaglia-Tsang ziggurat and Kinderman-Monahan-Leva ratio methods.
// The Ziggurat algorithm is the fastest available algorithm in most cases.

//double gsl_ran_gaussian_ziggurat(

// Normal random numbers, using the Ziggurat method
// This routine is based on the following article, with a couple of
// modifications which simplify the implementation.
//
//     George Marsaglia, Wai Wan Tsang
//     The Ziggurat Method for Generating Random Variables
//     Journal of Statistical Software, vol. 5 (2000), no. 8
//     http://www.jstatsoft.org/v05/i08/
//
// The modifications are:
//
// 1) use 128 steps instead of 256 to decrease the amount of static
// data necessary.
//
// 2) use an acceptance sampling from an exponential wedge
// exp(-R*(x-R/2)) for the tail of the base strip to simplify the
// implementation.  The area of exponential wedge is used in
// calculating 'v' and the coefficients in ziggurat table, so the
// coefficients differ slightly from those in the Marsaglia and Tsang
// paper.
//
// See also Leong et al, "A Comment on the Implementation of the
// Ziggurat Method", Journal of Statistical Software, vol 5 (2005), no 7.
//
// This function is based on the original GSL version, adapted to
// use R's system of RNGs. All credits to the original authors.
// Copyright (C) 2005  Jochen Voss.
// Copyright (C) 2020  Elias Haddad.
// License: GPL 3+
// GSL file: randist/gausszig.c
/// [[Rcpp::export]]
// Rcpp::NumericVector rnorm_zig(
//                             const gsl_rng * r
//                            ,const double sigma
//                           ){
//
//
//  /* position of right-most step */
//  #define PARAM_R 3.44428647676
//
//  unsigned long int i, j;
//  int sign;
//  double x, y;
//
//  const unsigned long int range = r->type->max - r->type->min;
//  const unsigned long int offset = r->type->min;
//
//  while (1)
//    {
//      if (range >= 0xFFFFFFFF)
//        {
//          unsigned long int k = gsl_rng_get(r) - offset;
//          i = (k & 0xFF);
//          j = (k >> 8) & 0xFFFFFF;
//        }
//      else if (range >= 0x00FFFFFF)
//        {
//          unsigned long int k1 = gsl_rng_get(r) - offset;
//          unsigned long int k2 = gsl_rng_get(r) - offset;
//          i = (k1 & 0xFF);
//          j = (k2 & 0x00FFFFFF);
//        }
//      else
//        {
//          i = gsl_rng_uniform_int (r, 256); /*  choose the step */
//          j = gsl_rng_uniform_int (r, 16777216);  /* sample from 2^24 */
//        }
//
//      sign = (i & 0x80) ? +1 : -1;
//      i &= 0x7f;
//
//      x = j * wtab[i];
//
//      if (j < ktab[i])
//        break;
//
//      if (i < 127)
//        {
//          double y0, y1, U1;
//          y0 = ytab[i];
//          y1 = ytab[i + 1];
//          U1 = gsl_rng_uniform (r);
//          y = y1 + (y0 - y1) * U1;
//        }
//      else
//        {
//          double U1, U2;
//          U1 = 1.0 - gsl_rng_uniform (r);
//          U2 = gsl_rng_uniform (r);
//          x = PARAM_R - log (U1) / PARAM_R;
//          y = exp (-PARAM_R * (x - 0.5 * PARAM_R)) * U2;
//        }
//
//      if (y < exp (-0.5 * x * x))
//        break;
//    }
//
//  return sign * sigma * x;
//}


//' Generates a Gamma-distributed sample
//'
//' This function returns a sample from a gamma-distributed random variable.
//' The method used to generate this sample is the Marsaglia-Tsang fast gamma
//' method. See more details below. The name was chosen so it doesn't clash
//' with R's native method.
//'
//' The gamma distribution is given by the function:
//' \deqn{f(x) = \frac{1}{\Gamma(b)a^b}x^{b-1}e^{-x/a}, x > 0}
//' where \eqn{b} is a shape parameter and \eqn{a} is a scale parameter.
//' The RNG is given by Marsaglia and Tsang, "A Simple Method for
//' generating gamma variables", ACM Transactions on Mathematical
//' Software, Vol 26, No 3 (2000), p363-372.
//' Available at:
//' https://doi.org/10.1145/358407.358414
//' This function is based on the original GSL version, adapted to
//' use R's system of RNGs. All credits to the original authors.
//' Implemented by J.D.Lamb@btinternet.com, minor modifications for GSL
//' by Brian Gough. Adapted to R by Elias Haddad.
//' Copyright (C) J.D.Lamb, Brian Gough.
//' Copyright (C) 2020-2021 Elias Haddad
//' License: GPL 3+
//' GSL file: randist/gamma.c
//'
//' @param n (int)
//' @param b (numeric) - shape parameter. Must be in the range \eqn{(0, \infty)}.
//' @param a (numeric) - scale parameter. Must be in the range \eqn{(0, \infty)}.
//' @return a numeric vector containing a random sample with above parameters.
//'
//' @examples
//' sample_gamma <- rgamma_c(1000, 1, 1)
//' @export
//' @md
// [[Rcpp::export]]
Rcpp::NumericVector rgamma_c(
  unsigned n
  ,double  b = 2
  ,double  a = 1/2
){

  // use R's RNG
  Rcpp::RNGScope scope;

  unsigned i=0;
  Rcpp::NumericVector sample(n);

  // assume a,b > 0
  if(b < 0 || a < 0){
    Rcpp::stop("The parameters a and b must be greater than zero.");
  }

  if( (b >= 0) & (b < 1) ){
    sample = Rcpp::runif(n); // identical to r
    return rgamma_c(n, 1.0 + b, a) * pow (sample, 1.0 / b);
  }else{

    double x, v;
    double d = b - 1.0 / 3.0;
    double c = (1.0 / 3.0) / sqrt (d);

    for(i=0; i<n; ++i){

      while(1){
        do{
          //x = gsl_ran_gaussian_ziggurat(r, 1.0);
          x = Rcpp::rnorm(1, 0, 1.0)[0];
          v = 1.0 + c * x;
        }while (v <= 0);

        v = v * v * v;
        sample[i] = Rcpp::runif(1)[0];

        if (sample[i] < 1 - 0.0331 * x * x * x * x)
          break;

        if (log (sample[i]) < 0.5 * x * x + d * (1 - v + log (v)))
          break;
      }

      sample[i] = a * d * v;
    }

    return sample;
    }
}

//' Generates a Laplace-distributed sample
//'
//' This function returns a sample from a Laplace-distributed random variable.
//'
//' The Laplace distribution is given by the two-sided exponential distribution
//' given by the function:
//' \deqn{ f(x;a,m) = \frac{1}{2a} e^{- \left| \frac{x-m}{a} \right|} }
//' The random sampling is done by inverse transform sampling.
//'
//' Copyright (C) 2020-2021 Elias Haddad
//' License: GPL 3+
//'
//' @param n (int) - the size of the sample.
//' @param m (numeric) - the location parameter.
//' @param a (numeric) - the scale parameter.
//' @return a numeric vector containing a random sample with above parameters.
//'
//' @examples
//' sample_gamma <- rlaplace(1000, 0, 1)
//' @export
//' @md
// [[Rcpp::export]]
Rcpp::NumericVector rlaplace(
  unsigned n
  ,const double m = 0
  ,const double a = 1
){

  // use R's RNG
  Rcpp::RNGScope scope;

  //Rcpp::NumericVector sample(n);
  Rcpp::NumericVector sample = Rcpp::runif(n);

  for(unsigned i=0; i<n; ++i){

    // this recovers the sign of the random draw
    double signal = sgn(sample[i]-0.5);

    // this calculates the inverse (quantile function)
    double value = m - signal*a*log(1 + signal*(1 - 2*sample[i]));

    sample[i] = value;
  }

  return sample;
}

//' Generates an Asymmetric Laplace-distributed sample
//'
//' This function returns a sample from an Asymetric Laplace distribution.
//'
//' The Asymmetric Laplace distribution is given by the two-sided exponential
//' distribution given by the function:
//' \deqn{f(x;a_l,a_r,m) =
//' \begin{cases}
//' \frac{1}{A} e^{-|\frac{x-m}{a_l}| }, & x < m
//' \frac{1}{A} e^{-|\frac{x-m}{a_r}| }, & x > m
//' \end{cases}}
//' with:
//' \deqn{A = a_l + a_r}
//' The random sampling is done by inverse transform sampling.
//'
//' Copyright (C) 2020-2021 Elias Haddad
//' License: GPL 3+
//'
//' @param n (int) - the size of the sample.
//' @param m (numeric) - the location parameter.
//' @param al,ar (numeric) - left and right scale parameters, respectively.
//' @return a numeric vector containing a random sample.
//'
//' @examples
//' sample_gamma <- ralaplace(1000)
//' @export
//' @md
// [[Rcpp::export]]
Rcpp::NumericVector ralaplace(
   unsigned n
  ,double   m  = 0
  ,double   al = 1
  ,double   ar = 1
){

  // use R's RNG
  Rcpp::RNGScope scope;

  //Rcpp::NumericVector sample(n);
  Rcpp::NumericVector sample = Rcpp::runif(n);
  double value;
  double signal;
  double A = al+ar;
  for(unsigned i=0; i<n; ++i){

    // this recovers the sign of the random draw
    signal = sgn(sample[i]- al/A);

    // this calculates the inverse (quantile function)
    //Rprintf("signal:%f", signal);
    if(signal < 0){
      value = m + al*log(A*sample[i]/al);
    }else{
      value = m - ar*log( A*(1 - sample[i])/ar );
    }
    sample[i] = value;
  }
  return sample;
}


//' Generates a random sample from a Exponential Power distribution
//'
//' This function returns a sample from a gamma-distributed random variable.
//'
//' The exponential power distribution (EP) is given by the function:
//' \deqn{ f(a,b) = \frac{1}{2a\Gamma(1+1/b)}e^{-|x/a|^b}, -\infty < x < \infty }.
//' where \eqn{b} is a shape parameter, \eqn{a} is a scale parameter and \eqn{\Gamma}
//' representes the gamma function. While not done here, this distribution can
//' be adapted to have non-zero location parameter.
//' The Exponential Power distribution is related to the gamma distribution by
//' the equation:
//' \deqn{E = a*G(1/b)^{1/b}}
//' where E and G are respectively EP and gamma random variables. This property
//' is used for cases where \eqn{b<1} and \eqn{b>4}. For \eqn{1 \leq b \leq 4}
//' rejection methods based on the Laplace and normal distributions are used,
//' which should be faster.
//' Technical details about this algorithm are available on:
//' P. R. Tadikamalla, "Random Sampling from the Exponential Power
//' Distribution", Journal of the American Statistical Association,
//' September 1980, Volume 75, Number 371, pages 683-686.
//' This function is based on the original GSL version, adapted to
//' use R's system of RNGs by Elias Haddad. All credits to the original authors.
//' Copyright (C) 1996, 1997, 1998, 1999, 2000, 2006, 2007 James Theiler, Brian Gough
//' Copyright (C) 2006 Giulio Bottazzi
//' Copyright (C) 2020-2021 Elias Haddad
//' License: GPL 3+
//' GSL file: randist/exppow.c
//' GSL function: gsl_ran_exppow
//'
//' @param n (int) - size of the sample.
//' @param m (numeric) - the location parameter.
//' @param a (numeric) - scale parameter.
//' @param b (numeric) - shape parameter.
//' @return a numeric vector containing a random sample with above parameters.
//'
//' @examples
//' sample_gamma <- rpower(1000)
//' @export
//' @md
// [[Rcpp::export]]
Rcpp::NumericVector rpower(
   unsigned n
  ,double   m = 0
  ,double   a = 1
  ,double   b = 2
){

  // use R's RNG
  Rcpp::RNGScope scope;

  Rcpp::NumericVector result(n);
  unsigned i=0;

  if (b < 1 || b > 4){

      // original implementation uses gsl_rng_uniform, which
      // allows for zero, while R's runif does not.
      // Does this impact the RNG in any meaningful way?
      Rcpp::NumericVector u = Rcpp::runif(n);
      Rcpp::NumericVector v = rgamma_c(n, 1 / b, 1.0);
      Rcpp::NumericVector z = a * pow (v, 1 / b);

      result = m + Rcpp::ifelse(u > 0.5, z, -z);
  }
  else if (b == 1){
      // Laplace distribution
    result =  rlaplace(n, m, a);
  }
  else if ( (b>1) && (b < 2) ){
      // Use laplace distribution for rejection method, from Tadikamalla

      double x, h, u;
      double B = pow (1 / b, 1 / b);

      for(i = 0; i<n; ++i){

        // rejection sampling
        do
          {
            x = rlaplace(1, 0, B)[0];
            u = Rcpp::runif(1)[0];
            h = -pow (fabs (x), b) + fabs (x) / B - 1 + (1 / b);
          }
        while (log (u) > h);

        // stores result in array
        result[i] = m + a*x;
      }
  }
  else if (b == 2){
    // Gaussian distribution
    result = Rcpp::rnorm(n, m, a/sqrt(2.0));
  }
  else if ( (b > 2) && (b < 4) ){
      // Use gaussian for rejection method, from Tadikamalla

      double x, h, u;
      double B = pow (1 / b, 1 / b);

      for(i = 0; i<n; ++i){

        // rejection sampling
        do
          {
            x = Rcpp::rnorm(1, 0, B)[0];
            u = Rcpp::runif(1)[0]; // this was the same RNG as R's
            h = -pow (fabs (x), b) + (x * x) / (2 * B * B) + (1 / b) - 0.5;
          }
        while (log (u) > h);

        // stores result in array
        result[i] = m + a*x;
      }
  }
  return result;
}


//' Produces a random sample from a Subbotin distribution
//'
//' Generate pseudo random-number from a Subbotin distribution using the
//' Tadikamalla method.
//'
//' The Subbotin distribution is given by the function:
//' \deqn{ f(x;a,b,m) = \frac{1}{A} e^{- \frac{1}{b} \left|\frac{x-m}{a}\right|^b} }
//' where \eqn{m} is a location parameter, \eqn{b} is a shape parameter, \eqn{a}
//' is a scale parameter and \eqn{\Gamma} representes the gamma function.
//' Since the Subbotin distribution is basically the exponential distribution
//' with scale parameter \eqn{a = ab^{1/b}} and \eqn{m=0}, we use the same
//' method of the exponential power RNG and add the location parameter.
//' Details can be found on the documentation of the \code{rpower} function.
//' Copyright (C) 2002-2014 Giulio Bottazzi
//' Copyright (C) 2020-2021 Elias Haddad
//' License: GPL 3+
//'
//' @param n (int) - the size of the sample.
//' @param m (numeric) - the location parameter.
//' @param a (numeric) - the scale parameter.
//' @param b (numeric) - the shape parameter.
//' @return a numeric vector containing a random sample.
//'
//' @examples
//' sample_gamma <- rsubbo(1000, 1, 1)
//' @export
//' @md
// [[Rcpp::export]]
Rcpp::NumericVector rsubbo(
   unsigned n
  ,double   m = 0
  ,double   a = 1
  ,double   b = 2
){

  // use R's RNG
  Rcpp::RNGScope scope;

 // calculate inverse distribution
  const double dtmp1 = a * pow(b, 1./b);
  Rcpp::NumericVector sample = m + rpower(n, 0, dtmp1, b);

  return sample;
}


//' Produces a random sample from a Asymmetric Power Exponential distribution
//'
//'
//' Generate pseudo random-number from an asymmetric power exponential distribution
//' using the Tadikamalla method.
//' This is the original version of Bottazzi (2004)
//'
//' The AEP distribution is expressed by the function:
//' \deqn{f(x;a_l,a_r,b_l,b_r,m) =
//' \begin{cases}
//' \frac{1}{A} e^{- \frac{1}{b_l} |\frac{x-m}{a_l}|^{b_l} }, & x < m
//' \frac{1}{A} e^{- \frac{1}{b_r} |\frac{x-m}{a_r}|^{b_r} }, & x > m
//' \end{cases} }
//' with:
//' \deqn{A = a_lb_l^{1/b_l}\Gamma(1+1/b_l) + a_rb_r^{1/b_r}\Gamma(1+1/b_r)}
//' where \eqn{m} is a location parameter, \eqn{b*} are shape parameters, \eqn{a*}
//' are scale parameters and \eqn{\Gamma} representes the gamma function.
//' By a suitably transformation, it is possible to use the EP distribution with
//' the Tadikamalla method to sample from this distribution. We basically take
//' the absolute values of the numbers sampled from the \code{rpower} function,
//' which is equivalent from sampling from a half Exponential Power distribution.
//' This values are then weighted by a constant expressed in the parameters.
//' More details are available on the package vignette and on the
//' function \code{rpower}.
//' Copyright (C) 2003-2014 Giulio Bottazzi
//' Copyright (C) 2020-2021 Elias Haddad
//' License: GPL 3+
//'
//' @param n (int) - size of the sample.
//' @param m (numeric) - location parameter.
//' @param bl,br (numeric) - shape parameters.
//' @param al,ar (numeric) - scale parameters.
//' @return a numeric vector containing a random sample.
//'
//' @examples
//' sample_gamma <- rasubbo(1000, 0, 0.5, 0.5,  1, 1)
//' @export
//' @md
// [[Rcpp::export]]
Rcpp::NumericVector rasubbo_orig(
   unsigned n
  ,double   m  = 0
  ,double   al = 1
  ,double   ar = 1
  ,double   bl = 2
  ,double   br = 2
){

  // use R's RNG
  Rcpp::RNGScope scope;

  // variables
  Rcpp::NumericVector sample(n);

  // generate variate
  const double Aleft  = al*pow(bl, 1./bl);
  const double Aright = ar*pow(br, 1./br);

  const double probleft = Aleft*gsl_sf_gamma(1./bl+1.)/
    (Aleft * gsl_sf_gamma(1./bl+1.) + Aright * gsl_sf_gamma(1./br+1.));

  sample = Rcpp::runif(n); // original implementation includes zero, while R's doesn't

  // calculates RNG
  for(unsigned i=0; i<n; ++i){

    if(sample[i] < probleft){
      sample[i] = m - fabs(rpower(1, 0, Aleft, bl)[0]);
    }
    else{
      sample[i] = m + fabs(rpower(1, 0, Aright, br)[0]);
    }
  }
  return sample;
}


//' Produces a random sample from a Asymmetric Power Exponential distribution
//'
//'
//' Generate pseudo random-number from an asymmetric power exponential distribution
//' using the Tadikamalla method.
//' This version improves on Bottazzi (2004) by making the mass of each
//' distribution to depend on the ratio between the \eqn{al} and the \eqn{ar}
//' parameters.
//'
//' The AEP distribution is expressed by the function:
//' \deqn{f(x;a_l,a_r,b_l,b_r,m) =
//' \begin{cases}
//' \frac{1}{A} e^{- \frac{1}{b_l} |\frac{x-m}{a_l}|^{b_l} }, & x < m
//' \frac{1}{A} e^{- \frac{1}{b_r} |\frac{x-m}{a_r}|^{b_r} }, & x > m
//' \end{cases} }
//' with:
//' \deqn{A = a_lb_l^{1/b_l}\Gamma(1+1/b_l) + a_rb_r^{1/b_r}\Gamma(1+1/b_r)}
//' where \eqn{m} is a location parameter, \eqn{b*} are shape parameters, \eqn{a*}
//' are scale parameters and \eqn{\Gamma} representes the gamma function.
//' By a suitably transformation, it is possible to use the EP distribution with
//' the Tadikamalla method to sample from this distribution. We basically take
//' the absolute values of the numbers sampled from the \code{rpower} function,
//' which is equivalent from sampling from a half Exponential Power distribution.
//' This values are then weighted by a constant expressed in the parameters.
//' More details are available on the package vignette and on the
//' function \code{rpower}.
//' Copyright (C) 2003-2014 Giulio Bottazzi
//' Copyright (C) 2020-2021 Elias Haddad
//' License: GPL 3+
//'
//' @param n (int) - size of the sample.
//' @param m (numeric) - location parameter.
//' @param bl,br (numeric) - shape parameters.
//' @param al,ar (numeric) - scale parameters.
//' @return a numeric vector containing a random sample.
//'
//' @examples
//' sample_gamma <- rasubbo(1000, 0, 0.5, 0.5,  1, 1)
//' @export
//' @md
// [[Rcpp::export]]
Rcpp::NumericVector rasubbo(
   unsigned n
  ,double   m  = 0
  ,double   al = 1
  ,double   ar = 1
  ,double   bl = 2
  ,double   br = 2
){

  // use R's RNG
  Rcpp::RNGScope scope;

  // variables
  Rcpp::NumericVector sample(n);
  double tmp;

  // generate variate
  // corretions here
  const double Aleft  = al*pow(bl, 1./bl)*gsl_sf_gamma(1./bl+1.);
  const double Aright = ar*pow(br, 1./br)*gsl_sf_gamma(1./br+1.);
  // this term creates the weighting
  const double Asum   = Aleft + Aright;

  const double probleft = Aleft*gsl_sf_gamma(1./bl+1.)/
    (Aleft * gsl_sf_gamma(1./bl+1.) + Aright * gsl_sf_gamma(1./br+1.));

  sample = Rcpp::runif(n); // original implementation includes zero, while R's doesn't

  // calculates RNG
  for(unsigned i=0; i<n; ++i){
    if(sample[i] < probleft){
      // here we use the weight to assess the mass in each side of the dist.
      tmp = bl*inv_inc_upper_gamma(1./bl, gsl_sf_gamma(1./bl)*Asum*sample[i]/Aleft);
      sample[i] = m - al*pow(tmp, 1./bl) ;
    }
    else{
      //sample[i] = m + fabs(rpower(1, 0, Aright, br)[0]);
      tmp = br*inv_inc_lower_gamma(1./br, gsl_sf_gamma(1./br)*Asum*(sample[i] - Aleft/Asum)/Aright);
      sample[i] = m + ar*pow(tmp, 1./br) ;
    }
  }
  return sample;
}
