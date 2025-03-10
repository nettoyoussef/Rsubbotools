/*
  laplafit (ver. 1.2.2) -- Fit an Asymmetric Laplace Distribution

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

/* Objective Function */
/*---------------- */

double lapla_nll(Rcpp::NumericVector data, const double m){

  unsigned size = data.size();
  unsigned utmp1;

  double sum=0.0;

  for(utmp1=0;utmp1<size;utmp1++)
    sum += fabs(m-data[utmp1]);

  return log(2*sum/size)+1.;

}

//' Fit a Laplace Distribution via maximum likelihood
//'
//' \code{laplafit} returns the parameters, standard errors. negative
//' log-likelihood and covariance matrix of the Laplace Distribution for a
//' sample. See details below.
//'
//' The Laplace distribution is a distribution controlled
//' by two parameters, with formula:
//' \deqn{f(x;a,m) = \frac{1}{2a} e^{- \left| \frac{x-m}{a} \right| }}
//' where \eqn{a} is a scale parameter, and \eqn{m} is a location parameter.
//' The estimations are produced by maximum likelihood, where analytical
//' formulas are available. Details on this method can be found on
//' the package vignette.
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
//' laplafit(sample_subbo)
//' @export
//' @md
// [[Rcpp::export]]
Rcpp::List laplafit(
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
  double m,a;

  // negative log-likelihood
  double fmin = 0;
  double tmp_fmin=0;

  // indices for pseudo-optimization
  long long int index, oldindex;
  long long int max=0, min=0;

  // index for general loops
  long long int i;

  /* store data */
  long long int size = data.size();/*the number of data*/

  /* sort data */
  /* --------- */
  sortRcpp(data);


  /* ---------------------- */

  /* store Fisher sqrt(variance-covariance)*/
  double sigma[3];


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
    fmin = lapla_nll(data, m);

    if(verb > 0){
      Rprintf("# index=%lld\n", index);
      Rprintf("#>>> Initial minimum: m=%e ll=%e\n", m, fmin);
    }

    // explore recursively the intervals on the right side
    do{

      oldindex = index;

      for(i=max+1; i <= max+interv_step && i < size-1; i++){

        tmp_fmin = lapla_nll(data, data[i]);

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

        tmp_fmin = lapla_nll(data, data[i]);

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
    Rprintf("#  intervals explored: %lld\n",max-min);
    Rprintf("\n");
  }

  /* decide the value of a  */
  a=0;
  for(i=0;i<size;i++){
    a+=fabs(m-data[i]);
  }
  a/=size;

  /* compute variance-covariance  */
  /* sigma[0] = std.err on m */
  /* sigma[1] = std.err on a */
  /* sigma[2] = corr.coeff m vs. a */
  sigma[0] = 2*a*a;
  sigma[1] = 1./sqrt(4.*a*a);
  sigma[2] = -1./sqrt(4.*a*a);

  // correlation matrix
  Rcpp::NumericMatrix I(2,2);
  I(0,0) = 1.;
  I(0,1) = sigma[2];
  I(1,0) = sigma[2];
  I(1,1) = 1.;


  // Name of parameters
  Rcpp::CharacterVector param_names = Rcpp::CharacterVector::create("m", "a");

  Rcpp::NumericVector par =
    Rcpp::NumericVector::create(m ,a);

  Rcpp::NumericVector std_error =
    Rcpp::NumericVector::create(sigma[0]/sqrt(size), sigma[1]/sqrt(size));

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
