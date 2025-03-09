/*
  alaplafit (ver. 1.2.1) -- Fit an Asymmetric Laplace Distribution

  Copyright (C) 2004-2014 Giulio Bottazzi
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

/*----------------- */
/* Object Function */
/*---------------- */

double alapla_nll(Rcpp::NumericVector data, const double m){

  unsigned size = data.size();
  unsigned utmp1;

  double sumL=0.0;
  double sumR=0.0;


/*   fprintf(Fmessages,"#objf bl=%.3e br=%.3e al=%.3e ar=%.3e m=%.3e\n", */
/*        x[0],x[1],x[2],x[3],x[4]); */

  for(utmp1=0;utmp1<size;utmp1++){
    if(data[utmp1]>m) break;
    sumL+=m-data[utmp1];
  }
  for(;utmp1<size;utmp1++){
    sumR+=data[utmp1]-m;
  }

  sumL = sqrt(sumL/size);
  sumR = sqrt(sumR/size);

  return 2*log(sumL+sumR)+1.;

}

//' Fit an Asymmetric Laplace Distribution via maximum likelihood
//'
//' \code{alaplafit} returns the parameters, standard errors. negative
//' log-likelihood and covariance matrix of the Asymmetric Laplace Distribution
//' for a sample. See details below.
//'
//' The Asymmetric Laplace distribution is a distribution controlled
//' by three parameters, with formula:
//' \deqn{f(x;a_l,a_r,m) =
//' \begin{cases}
//' \frac{1}{A} e^{-|\frac{x-m}{a_l}| }, & x < m \\
//' \frac{1}{A} e^{-|\frac{x-m}{a_r}| }, & x > m
//' \end{cases}}
//' with:
//' \deqn{A = a_l + a_r}
//' where \eqn{a*} are scale parameters, and \eqn{m} is a location parameter.
//' It is basically derived from the Asymmetric Exponential Power distribution
//' by setting \eqn{b_l = b_r = b}.
//' The estimations are produced by maximum likelihood, where
//' analytical formulas are available for the \eqn{a*} parameters.
//' The \eqn{m} parameter is found by an iterative method, using the median as
//' the initial guess. The method explore intervals around the last minimum
//' found, similar to the \code{subboafit} routine.
//' Details on this method can be found on the package vignette.
//'
//' @param data (NumericVector) - the sample used to fit the distribution.
//' @param verb (int) - the level of verbosity. Select one of:
//' * 0  just the final result
//' * 1  details of optim. routine
//' @param interv_step  int - the number of intervals to be explored after
//' the last minimum was found in the interval optimization. Default is 10.
//' @param provided_m_ NumericVector - if NULL, the m parameter is estimated
//' by the routine. If numeric, the estimation fixes m to the given value.
//' @return a list containing the following items:
//' * "dt" - dataset containing parameters estimations and standard deviations.
//' * "log-likelihood" - negative log-likelihood value.
//' * "matrix" - the covariance matrix for the parameters.
//'
//' @examples
//' sample_subbo <- rpower(1000, 1, 1)
//' alaplafit(sample_subbo)
//' @export
//' @md
// [[Rcpp::export]]
Rcpp::List alaplafit(
                     Rcpp::NumericVector data
                     ,int verb = 0
                     ,unsigned interv_step = 10
                     ,Rcpp::Nullable<Rcpp::NumericVector> provided_m_ = R_NilValue
                     ){


  /* initial values */
  /* -------------- */

  /* store possibly provided values for parameters */
  unsigned is_m_provided = 0;

  // initial value for m parameter
  if(provided_m_.isNotNull()){

    is_m_provided = 1;

    }

  /* store guess */
  double m,al,ar;

  /* store loglike */
  double fmin = 0;
  double tmp_fmin;

  // indices for pseudo-optimization
  size_t index, oldindex, max, min;

  // index for general loops
  size_t i;

  /* store data */
  size_t size = data.size();/*the number of observations*/

  // partial sums for the negative likelihood
  double sl=0,sr=0;

  /* sort data */
  /* --------- */
  sortRcpp(data);

  /* store Fisher sqrt(variance-covariance)*/
  double sigma[5];


  /* decide the value of m  */
  if (is_m_provided){
    // casting the null value
    // necessary according to
    // https://gallery.rcpp.org/articles/optional-null-function-arguments/
    Rcpp::NumericVector provided_m(provided_m_);
    m = provided_m[0];
  }
  else{

    index = calculate_index(size);

    max = min = index;
    m = median(data, size);
    fmin = alapla_nll(data, m);

    if(verb > 0){
      Rprintf("# index=%ld\n", index);
      Rprintf("#>>> Initial minimum: m=%e ll=%e\n", m, fmin);
    }

    // explore recursively the intervals on the right side
    do{

      oldindex = index;

      for(i=max+1; i <= max+interv_step && i < size-1; i++){

        tmp_fmin = alapla_nll(data, data[i]);

        if(tmp_fmin < fmin){/* found new minimum */
          m    = data[i];
          fmin = tmp_fmin;
          index = i;

          // print results - we use arrows to
          // differentiate the global minimum
          if(verb > 0){
            Rprintf("#>>> [%+.3e:%+.3e] m=%e ll=%e\n"
                    , data[i], data[i+1], m, tmp_fmin);
          }
        }else{
          // print results
          if(verb > 0){
            Rprintf("#    [%+.3e:%+.3e] ll=%e\n"
                    , data[i], data[i+1], tmp_fmin);
          }
        }
      }
      // counts the number of intervals explored
      // and adjusts max for the next iteration, if it occurs
      max = i-1;

    }while(oldindex != index);

    //  // explore recursively the intervals on the left side
    do{

      oldindex = index;

      for(i=min-1; i >= min-interv_step && i >= 0; i--){

        tmp_fmin = alapla_nll(data, data[i]);

        if(tmp_fmin < fmin){/* found new minimum */
          m     = data[i];
          fmin  = tmp_fmin;
          index = i;

          // print results - we use arrows to
          // differentiate the global minimum
          if(verb > 0){
            Rprintf("#>>> [%+.3e:%+.3e] m=%e ll=%e\n"
                    , data[i], data[i+1], m, tmp_fmin);
          }
        }else{
          // print results
          if(verb > 0){
            Rprintf("#    [%+.3e:%+.3e] ll=%e\n"
                    , data[i], data[i+1], tmp_fmin);
          }
        }

      }
      // counts the number of intervals explored
      // and adjusts min for the next iteration, if it occurs
      min = i+1;

    }while(oldindex != index);
  }

  if(verb > 0){
    Rprintf("Results of interval optimization: \n");
    Rprintf("#>>> m=%e ll=%e\n", m, tmp_fmin);
    Rprintf("\n");
    Rprintf("#  intervals explored: %ld\n",max-min);
    Rprintf("\n");
  }



  /* decide the value of al and ar  */
  for(i=0;i<size;i++){
    if(data[i]>m) break;
    sl+=m-data[i];
  }
  sl/=size;

  for(;i<size;i++){
    sr+=data[i]-m;
  }
  sr/=size;

  al=sl+sqrt(sl*sr);
  ar=sr+sqrt(sl*sr);

  /* compute variance-covariance  */
  /* sigma[0] = std.err on m */
  /* sigma[1] = std.err on al */
  /* sigma[2] = std.err on ar */
  /* sigma[3] = corr.coeff m vs. al */
  /* sigma[4] = corr.coeff m vs. ar */
  sigma[0] = 2*al*ar;
  sigma[1] = al*(al+ar);
  sigma[2] = ar*(al+ar);
  sigma[3] = 1./sqrt(2.*ar*(al+ar));
  sigma[4] = -1./sqrt(2.*al*(al+ar));

  // correlation matrix
  Rcpp::NumericMatrix I(3,3);
  I(0,0) = I(1,1) = I(2,2) = 1.;
  I(0,1) = I(1,0) = sigma[3];
  I(0,2) = I(2,0) = sigma[4];
  I(1,2) = I(2,1) = 0.;

  // Name of parameters
  Rcpp::CharacterVector param_names =
    Rcpp::CharacterVector::create("m", "al", "ar");

  Rcpp::NumericVector par =
    Rcpp::NumericVector::create(m ,al, ar);

  Rcpp::NumericVector std_error =
    Rcpp::NumericVector::create(sigma[0]/sqrt(size)
                                ,sigma[1]/sqrt(size)
                                ,sigma[2]/sqrt(size)
                                );

  // main dataframe with coefficients
       Rcpp::List dt =
       Rcpp::DataFrame::create( Rcpp::Named("param")      = param_names
                                ,Rcpp::Named("coef")      = par
                                ,Rcpp::Named("std_error") = std_error
                                );

       Rcpp::List ans = Rcpp::List::create(
                                      Rcpp::Named("dt")              = dt
                                      ,Rcpp::Named("log-likelihood") = fmin
                                      ,Rcpp::Named("matrix")         = I
                                      );

  return ans;
}
