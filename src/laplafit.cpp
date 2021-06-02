/*
  laplafit (ver. 1.2.2) -- Fit an asymmetric exponential density

  Copyright (C) 2004-2014 Giulio Bottazzi

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


/*
  Verbosity levels:
  0 just the final ouput
  1 the results of intermediate steps
  2 internal information on intermediate steps
  3 gory details on intermediate steps
*/


#include "common.h"


/* Global Variables */
/* ---------------- */



/*------------------*/


/* Output Functions */
/*----------------- */

void lapla_printcumul(Rcpp::NumericVector data, const double m, const double a){

  int size = data.size();
  int i;
  double dtmp1;

  for(i=0;i<size;i++){
    if(data[i]>m){
      dtmp1=1-exp( (m-data[i])/a )/2;
    }
    else{
      dtmp1=exp( (data[i]-m)/a )/2;
    }
    Rprintf("%e %e\n",data[i],dtmp1);
  }

}

void lapla_printdensity(Rcpp::NumericVector data, const double m, const double a){

  int size = data.size();
  int i;

  for(i=0;i<size;i++){
    double dtmp1=data[i];
    Rprintf("%e ",dtmp1);
    dtmp1=dtmp1-m;
    if(dtmp1>=0){
      Rprintf("%e\n",exp(-dtmp1/a)/(2*a));
    }
    else{
      Rprintf("%e\n",exp( dtmp1/a)/(2*a));
    }
  }

}
/*----------------- */


/* Objective Function */
/*---------------- */

double lapla_nll(Rcpp::NumericVector data, const double m){

  int size = data.size();
  unsigned utmp1;

  double sum=0.0;

  for(utmp1=0;utmp1<size;utmp1++)
    sum += fabs(m-data[utmp1]);

  return log(2*sum/size)+1.;

}





// Fit Laplace density. Read data from files or from standard input.  \n");
// With output option '4' prints the log-likelihood as a function of m   \n");
// Usage: %s [options] [files]\n\n",argv[0]);
// input.                                                             \n\n");
//  Options:                                                            \n");
//  -O  output type (default 0)                \n");
//       0  parameters m a and negative log-likelihood\n");
//       1  the estimated distribution function computed on the provided points     \n");
//       2  the estimated density function computed on the provided points \n");
//       3  parameters m, a and their standard errors and correlations \n");
//       4  log-likelihood profile \n");
//  -m  the mode is not estimated but is set to the value provided\n");
//  -x  initial value of m or plot range if output type is '4' (default 0)\n");
//  -n  number of plotted points if output type is '4' (default 10)\n");
// Examples:\n");
//  'laplafit -m 1 <file'  estimate a with m=1\n");
//  'laplafit -O 4 -x -1,1 <file'  print the log-likelihood on a regular grid\n");
// from -1 to 1. The grid has 10 points. \n");
// [[Rcpp::export]]
Rcpp::List laplafit(
                    Rcpp::NumericVector data
                    ,int verb = 0
                    ,int method = 7
                    ,int output = 0
                    ,Rcpp::Nullable<Rcpp::NumericVector> provided_m_ = R_NilValue
                    ,unsigned interv_step = 10
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
  double tmp_fmin;

  // indices for pseudo-optimization
  size_t index, oldindex, max, min;

  // index for general loops
  size_t i;

  /* store data */
  size_t size = data.size();/*the number of data*/

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

    if(verb > 1){
      Rprintf("# index=%d\n", index);
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
          if(verb > 1){
            Rprintf("#>>> [%+.3e:%+.3e] m=%e ll=%e\n"
                    , data[i], data[i+1], m, tmp_fmin);
          }
        }else{
          // print results
          if(verb > 1){
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
          if(verb > 1){
            Rprintf("#>>> [%+.3e:%+.3e] m=%e ll=%e\n"
                    , data[i], data[i+1], m, tmp_fmin);
          }
        }else{
          // print results
          if(verb > 1){
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

  if(verb > 1){
    Rprintf("Results of interval optimization: \n");
    Rprintf("#>>> m=%e ll=%e\n", m, tmp_fmin);
    Rprintf("\n");
    Rprintf("#  intervals explored: %d\n",max-min);
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
