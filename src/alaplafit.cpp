/*
  alaplafit (ver. 1.2.1) -- Fit an asymmetric exponential density

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




/* Output Functions */
/*----------------- */

void alapla_printcumul(Rcpp::NumericVector data, const double m, const double al, const double ar){

  int size = data.size();
  int i;
  double dtmp1;

  const double Asum=al+ar;

  for(i=0;i<size;i++){
    if(data[i]>m){
      dtmp1=(data[i]-m)/ar;
      dtmp1=1-exp(-dtmp1)*ar/Asum;
     }
    else{
      dtmp1=(m-data[i])/al;
      dtmp1=exp(-dtmp1)*al/Asum;
    }
    Rprintf("%e %e\n",data[i],dtmp1);
  }

}

void alapla_printdensity(Rcpp::NumericVector data, const double m, const double al, const double ar){

  int size = data.size();
  int i;

  const double norm=al+ar;

  for(i=0;i<size;i++){
    double dtmp1=data[i];
    Rprintf("%e ",dtmp1);
    dtmp1=dtmp1-m;
    if(dtmp1>=0){
      Rprintf("%e\n",exp(-dtmp1/ar)/norm);
    }
    else{
      Rprintf("%e\n",exp( dtmp1/al)/norm);
    }
  }

}
/*----------------- */


/* Object Function */
/*---------------- */

double alapla_nll(Rcpp::NumericVector data, const double m){

  int size = data.size();
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


// Fit asymmetric Laplace density. Read data from files or from standard input  \n");
// With output option '4' prints the log-likelihood as a function of m   \n");
// Usage: %s [options] [files]\n\n",argv[0]);
// input.                                                             \n\n");
//  Options:                                                            \n");
//  -O  output type (default 0)                \n");
//       0  parameters m al ar and negative log-likelihood\n");
//       1  the estimated distribution function computed on the provided points     \n");
//       2  the estimated density function computed on the provided points \n");
// 3  parameters m,al,ar their standard errors and correlations \n");
//       4  log-likelihood profile \n");
//  -m  the mode is not estimated but is set to the value provided\n");
//  -x  initial value of m or plot range if output type is '4' (default 0)\n");
//  -n  number of plotted points if output type is '4' (default 10)\n");
// Examples:\n");
//  'alaplafit -m 1 <file'  estimate al and ar with m=1\n");
//  'alaplafit -O 4 -x -1,1 <file'  print the log-likelihood on a regular grid\n");
// from -1 to 1. The grid has 10 points. \n");
// [[Rcpp::export]]
Rcpp::List alaplafit(
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

    if(verb > 1){
      Rprintf("# index=%d\n", index);
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

        tmp_fmin = alapla_nll(data, data[i]);

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
