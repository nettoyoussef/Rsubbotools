/*
  subboafit (ver. 1.2.1) -- Fit an asymmetric power exp. density via likelihood maximization

  Copyright (C) 2007 Giulio Bottazzi
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

void subboa_printcumul(Rcpp::NumericVector data, double param[]){

  int size = data.size();
  int i;
  double dtmp1;
  const double bl=param[0];
  const double br=param[1];
  const double al=param[2];
  const double ar=param[3];
  const double m=param[4];

  const double Al=al*gsl_sf_gamma(1./bl+1.)*pow(bl,1./bl);
  const double Ar=ar*gsl_sf_gamma(1./br+1.)*pow(br,1./br);
  const double Asum=Al+Ar;

  for(i=0;i<size;i++){
    if(data[i]>m){
      dtmp1=pow((data[i]-m)/ar,br)/br;
      dtmp1=(Al+Ar*gsl_sf_gamma_inc_P(1./br,dtmp1))/Asum;
     }
    else{
      dtmp1=pow((m-data[i])/al,bl)/bl;
      dtmp1=Al*gsl_sf_gamma_inc_Q(1./bl,dtmp1)/Asum;
    }
    Rprintf("%e %e\n",data[i],dtmp1);
  }

}

void subboa_printdensity(Rcpp::NumericVector data, double param[]){

  int size = data.size();
  int i;
  const double bl=param[0];
  const double br=param[1];
  const double al=param[2];
  const double ar=param[3];
  const double m=param[4];

  const double norm=al*pow(bl,1/bl)*gsl_sf_gamma(1./bl+1.)
    +ar*pow(br,1/br)*gsl_sf_gamma(1./br+1);

  for(i=0;i<size;i++){
    double dtmp1=data[i];
    Rprintf("%e ",dtmp1);
    dtmp1-=m;
    if(dtmp1>=0){
      Rprintf("%e\n",exp(-pow(dtmp1/ar,br)/br)/norm);
    }
    else{
      Rprintf("%e\n",exp(-pow(-dtmp1/al,bl)/bl)/norm);
    }
  }

}
/*----------------- */



/*
   par   array of paramaters
   N     number of observations
   dim   dimension of the matrix: 2 m known; 3 m unknown
   I     the variance-covarance matrxi
*/
RcppGSL::Matrix subboa_varcovar(const Rcpp::NumericVector par, const size_t N, const size_t dim){

  const double bl = par[0];
  const double br = par[1];
  const double al = par[2];
  const double ar = par[3];
  //const double m  = par[4]; // declared but not used in the original

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

  size_t i,j;

  // declare matrix
  RcppGSL::Matrix I(dim, dim);

  // needed for inversion
  RcppGSL::Matrix J(dim,dim);
  gsl_permutation *P = gsl_permutation_alloc (dim);
  int signum;


  // define the matrix, the order is bl,br,al,ar,m

  // bl - bl
  I(0,0) = al*(dB0ldx2-al*dB0ldx*dB0ldx/A +
               B2l/bl -2*B1l/(bl*bl)+2*B0l/pow(bl,3))/A ;

  // bl - br
  // apparently I can't do this in one line
  //I(0,1) = I(1,0) = -al*ar*dB0ldx*dB0rdx/(A*A);
  // in fact, there is a bug with RcppGSL
  // I can't use the format I(x,t) = I(y,k)
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

  // al - al
  I(2,2) = -B0l*B0l/(A*A)+(bl+1)*B0l/(al*A);

  // al - ar
  I(2,3) = -B0l*B0r/(A*A);
  I(3,2) = gsl_matrix_get(I,2,3);

  // ar - ar
  I(3,3) = -B0r*B0r/(A*A)+(br+1)*B0r/(ar*A);


  if(dim == 5){

    const double dt1l = gsl_sf_gamma (2.-1/bl);
    const double dt1r = gsl_sf_gamma (2.-1/br);
    const double dt2l = pow(bl,1.-1./bl);
    const double dt2r = pow(br,1.-1./br);

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

    // m - m
    I(4,4) = ( dt1l*dt2l/al+dt1r*dt2r/ar)/A ;

  }

  // print matrix (debugging purposes)
  //for(i = 0;i<dim;i++){
  //  for(j = 0;j<dim;j++){
  //    Rprintf("%10.10f ", gsl_matrix_get(I,i,j));
  //  }
  //  Rprintf("\n");
  //}
  //Rprintf("-------------\n");

  // invert I; in J store temporary LU decomp.
  gsl_matrix_memcpy (J,I);
  gsl_linalg_LU_decomp (J,P,&signum);
  gsl_linalg_LU_invert (J,P,I);

  // free allocated memory
  gsl_permutation_free(P);

  // print matrix
  //for(i = 0;i<dim;i++){
  //  for(j = 0;j<dim;j++){
  //    Rprintf("%10.10f ", gsl_matrix_get(I,i,j));
  //  }
  //  Rprintf("\n");
  //}
  //Rprintf("-------------\n");

  // set the var-covar matrix
  for(i=0;i<dim;i++){
    for(j=0;j<i;j++){
      I(i,j) =
        gsl_matrix_get(I,i,j)/sqrt(gsl_matrix_get(I,i,i)*gsl_matrix_get(I,j,j));
    }

  }

  // print matrix
  //for(i = 0;i<dim;i++){
  //  for(j = 0;j<dim;j++){
  //    Rprintf("%10.10f ", gsl_matrix_get(I,i,j));
  //  }
  //  Rprintf("\n");
  //}
  //Rprintf("-------------\n");

  for(i=0;i<dim;i++){
    for(j=i;j<dim;j++){
      I(i,j) = gsl_matrix_get(I,i,j)/N;
    }

  }

  // print matrix
  //for(i = 0;i<dim;i++){
  //  for(j = 0;j<dim;j++){
  //    Rprintf("%10.10f ", gsl_matrix_get(I,i,j));
  //  }
  //  Rprintf("\n");
  //}
  //Rprintf("-------------\n");

 return I;

}



/* Objective Function */
/*---------------- */

void subboa_objf(Rcpp::NumericVector data, const size_t n, Rcpp::NumericVector x, void *params, double *f){

  // data size
  int size = data.size();

  // temporary variables
  double dtmp1;
  unsigned utmp1;

  double sumL=0.0;
  double sumR=0.0;

  // deferencing variables to
  // make code more clear
  const double bl=x[0];
  const double br=x[1];
  const double al=x[2];
  const double ar=x[3];
  const double mu=x[4];

/*   fprintf(Fmessages,"#objf bl=%.3e br=%.3e al=%.3e ar=%.3e m=%.3e\n", */
/*        x[0],x[1],x[2],x[3],x[4]); */

  for(utmp1=0;utmp1<size;utmp1++){
    if(data[utmp1]<mu){
      sumL+=pow(mu-data[utmp1],bl);
    }
    else if (data[utmp1]>mu){
      sumR+=pow(data[utmp1]-mu,br);
    }
  }

/*   for(utmp1=0;utmp1<size;utmp1++){ */
/*     if(data[utmp1]>mu) break;     */
/*     sumL+=pow(mu-data[utmp1],bl); */
/*   } */
/*   for(;utmp1<size;utmp1++){ */
/*     sumR+=pow(data[utmp1]-mu,br); */
/*   } */

  sumL /= size;
  sumR /= size;

  dtmp1 = al*pow(bl,1./bl)*gsl_sf_gamma(1./bl+1.)
    +ar*pow(br,1./br)*gsl_sf_gamma(1./br+1.);

  *f = log(dtmp1)
    +sumL/(pow(al,bl)*bl)+sumR/(pow(ar,br)*br);

}



void subboa_objdf(Rcpp::NumericVector data, const size_t n, Rcpp::NumericVector x, void *params, Rcpp::NumericVector df){


  // data size
  int size = data.size();

  // temp variables
  double dtmp1;

  // deferencing variables to
  // make code more clear
  const double bl=x[0];
  const double br=x[1];
  const double al=x[2];
  const double ar=x[3];
  const double mu=x[4];

  unsigned utmp1;
  double sumL=0.0;
  double sumR=0.0;
  double sumL1=0.0;
  double sumR1=0.0;
  double sumLl=0.0;
  double sumRl=0.0;
/*   fprintf(Fmessages,"#objdf bl=%.3e br=%.3e al=%.3e ar=%.3e m=%.3e\n", */
/*        x[0],x[1],x[2],x[3],x[4]); */



  for(utmp1=0;utmp1<size;utmp1++){
    if(data[utmp1]<mu){
      sumL1+= (dtmp1 = pow(mu-data[utmp1],bl-1.) );
      dtmp1*= mu-data[utmp1];
      sumL+= dtmp1;
      sumLl+=log(mu-data[utmp1])*dtmp1;
    }
    else if (data[utmp1]>mu){
      sumR1+=(dtmp1 = pow(data[utmp1]-mu,br-1.));
      dtmp1*= data[utmp1]-mu;
      sumR+=dtmp1;
      sumRl+=log(data[utmp1]-mu)*dtmp1;
    }
  }

  sumL /= size;
  sumR /= size;
  sumL1 /= size;
  sumR1 /= size;
  sumLl /= size;
  sumRl /= size;

  dtmp1 = al*pow(bl,1./bl)*gsl_sf_gamma(1./bl+1.)
    +ar*pow(br,1./br)*gsl_sf_gamma(1./br+1.);

  df[0] = (1.-log(bl)-gsl_sf_psi(1./bl+1.))*gsl_sf_gamma(1./bl+1.)*pow(bl,1./bl-2.)*al/dtmp1
    -(1./(bl*bl) + log(al)/bl)*sumL/pow(al,bl) + sumLl/(pow(al,bl)*bl);

  df[1] = (1.-log(br)-gsl_sf_psi(1./br+1.))*gsl_sf_gamma(1./br+1.)*pow(br,1./br-2.)*ar/dtmp1
    -(1./(br*br) + log(ar)/br)*sumR/pow(ar,br) + sumRl/(pow(ar,br)*br);

  df[2] = pow(bl,1./bl)*gsl_sf_gamma(1./bl+1.)/dtmp1
    - sumL/pow(al,bl+1.);

  df[3] = pow(br,1./br)*gsl_sf_gamma(1./br+1.)/dtmp1
    - sumR/pow(ar,br+1.);

  df[4] = sumL1/pow(al,bl) - sumR1/pow(ar,br);


}


void subboa_objfdf(Rcpp::NumericVector data, const size_t n, Rcpp::NumericVector x, void *params, double *f, Rcpp::NumericVector df){

  // data size
  int size = data.size();

  // temp variables
  double dtmp1;

  // deferencing variables to
  // make code more clear
  const double bl=x[0];
  const double br=x[1];
  const double al=x[2];
  const double ar=x[3];
  const double mu=x[4];

  unsigned utmp1;
  double sumL=0.0;
  double sumR=0.0;
  double sumL1=0.0;
  double sumR1=0.0;
  double sumLl=0.0;
  double sumRl=0.0;

/*   fprintf(Fmessages,"#objfdf bl=%.3e br=%.3e al=%.3e ar=%.3e m=%.3e\n", */
/*        x[0],x[1],x[2],x[3],x[4]); */

  for(utmp1=0;utmp1<size;utmp1++){
    if(data[utmp1]<mu){
      sumL1+= (dtmp1 = pow(mu-data[utmp1],bl-1.) );
      dtmp1*= mu-data[utmp1];
      sumL+= dtmp1;
      sumLl+=log(mu-data[utmp1])*dtmp1;
    }
    else if (data[utmp1]>mu){
      sumR1+=(dtmp1 = pow(data[utmp1]-mu,br-1.));
      dtmp1*= data[utmp1]-mu;
      sumR+=dtmp1;
      sumRl+=log(data[utmp1]-mu)*dtmp1;
    }
  }

  sumL /= size;
  sumR /= size;
  sumL1 /= size;
  sumR1 /= size;
  sumLl /= size;
  sumRl /= size;

  dtmp1 = al*pow(bl,1./bl)*gsl_sf_gamma(1./bl+1.)
    +ar*pow(br,1./br)*gsl_sf_gamma(1./br+1.);

  *f = log(dtmp1)
    +sumL/(pow(al,bl)*bl)+sumR/(pow(ar,br)*br);

  df[0] = (1.-log(bl)-gsl_sf_psi(1./bl+1.))*gsl_sf_gamma(1./bl+1.)*pow(bl,1./bl-2.)*al/dtmp1
    -(1./(bl*bl) + log(al)/bl)*sumL/pow(al,bl) + sumLl/(pow(al,bl)*bl);

  df[1] = (1.-log(br)-gsl_sf_psi(1./br+1.))*gsl_sf_gamma(1./br+1.)*pow(br,1./br-2.)*ar/dtmp1
    -(1./(br*br) + log(ar)/br)*sumR/pow(ar,br) + sumRl/(pow(ar,br)*br);

  df[2] = pow(bl,1./bl)*gsl_sf_gamma(1./bl+1.)/dtmp1
    - sumL/pow(al,bl+1.);

  df[3] = pow(br,1./br)*gsl_sf_gamma(1./br+1.)/dtmp1
    - sumR/pow(ar,br+1.);

  df[4] = sumL1/pow(al,bl) - sumR1/pow(ar,br);

}
/*---------------- */


//' Fit an Asymmetric Power Exponential density via maximum likelihood
//'
//' \code{subboafit} returns the parameters, standard errors. negative
//' log-likelihood and covariance matrix of the asymmetric power exponential for
//' a sample. The process can execute two steps, dependending on the level of
//' accuracy required. See details below.
//'
//' The AEP is a exponential power distribution controlled
//' by five parameters, with formula:
//' \deqn{ f(x;a_l,a_r,b_l,b_r,m) =
//' \begin{cases}
//' \frac{1}{A} e^{- \frac{1}{b_l} |\frac{x-m}{a_l}|^{b_l} }, & x < m \\
//' \frac{1}{A} e^{- \frac{1}{b_r} |\frac{x-m}{a_r}|^{b_r} }, & x > m
//' \end{cases} }
//' with:
//' \deqn{A = a_lb_l^{1/b_l}\Gamma(1+1/b_l) + a_rb_r^{1/b_r}\Gamma(1+1/b_r)}
//' where \eqn{l} and \eqn{r} represent left and right tails, \eqn{a*} are
//' scale parameters, \eqn{b*} control the tails (lower values represent
//' fatter tails), and \eqn{m} is a location parameter. Due to its lack of
//' simmetry, and differently from the Subbotin, there is no simple equations
//' available to use the method of moments, so we start directly by minimizing
//' the negative log-likelihood. This global optimization is executed without
//' restricting any parameters. If required (default), after the global
//' optimization is finished, the method proceeds to iterate over the intervals
//' between several two observations, iterating the same algorithm of the
//' global optimization. The last method happens because of the lack of
//' smoothness on the \eqn{m} parameter, and intervals must be used since the
//' likelihood function doesn't have a derivative whenever \eqn{m} equals a
//' sample observation. Due to the cost, these iterations are capped at most
//' \emph{interv_step} (default 10) from the last minimum observed.
//' Details on this method are available on the package vignette.
//'
//' @param data (NumericVector) - the sample used to fit the distribution.
//' @param verb (int) - the level of verbosity. Select one of:
//' * 0  just the final result
//' * 1  headings and summary table
//' * 2  intermediate steps results
//' * 3  intermediate steps internals
//' * 4+  details of optim. routine
//' @param method int - the steps that should be used to estimate the
//' parameters.
//' * 0 no optimization perform - just return the log-likelihood from initial guess.
//' * 1 global optimization not considering lack of smoothness in m
//' * 2 interval optimization taking non-smoothness in m into consideration
//' @param interv_step  int - the number of intervals to be explored after
//' the last minimum was found in the interval optimization. Default is 10.
//' @param provided_m_ NumericVector - if NULL, the m parameter is estimated
//' by the routine. If numeric, the estimation fixes m to the given value.
//' @param par NumericVector - vector containing the initial guess for
//' parameters bl, br, al, ar and m, respectively. Default values of are
//' c(2, 2, 1, 1, 0).
//' @param g_opt_par NumericVector - vector containing the global optimization
//' parameters.
//' The optimization parameters are:
//' * step  - (num) initial step size of the searching algorithm.
//' * tol   - (num) line search tolerance.
//' * iter  - (int) maximum number of iterations.
//' * eps   - (num) gradient tolerance. The stopping criteria is \eqn{||\text{gradient}||<\text{eps}}.
//' * msize - (num) simplex max size. stopping criteria given by \eqn{||\text{max edge}||<\text{msize}}
//' * algo  - (int) algorithm. the optimization method used:
//'   * 0 Fletcher-Reeves
//'   * 1 Polak-Ribiere
//'   * 2 Broyden-Fletcher-Goldfarb-Shanno
//'   * 3 Steepest descent
//'   * 4 Nelder-Mead simplex
//'   * 5 Broyden-Fletcher-Goldfarb-Shanno ver.2
//'
//' Details for each algorithm are available on the [GSL Manual](https://www.gnu.org/software/gsl/doc/html/multimin.html).
//' Default values are c(.1, 1e-2, 100, 1e-3, 1e-5, 2).
//' @param itv_opt_par NumericVector - interval optimization parameters. Fields
//' are the same as the ones for the global optimization. Default values
//' are c(.01, 1e-3, 200, 1e-3, 1e-5, 5).
//' @return a list containing the following items:
//' * "dt" - dataset containing parameters estimations and standard deviations.
//' * "log-likelihood" - negative log-likelihood value.
//' * "matrix" - the covariance matrix for the parameters.
//'
//' @examples
//' sample_subbo <- rpower(1000, 1, 2)
//' subboafit(sample_subbo)
//' @export
//' @md
// [[Rcpp::export]]
Rcpp::List subboafit(
                     Rcpp::NumericVector data
                     ,int verb = 0
                     ,int method = 6
                     ,int interv_step = 10
                     ,Rcpp::Nullable<Rcpp::NumericVector> provided_m_ = R_NilValue
                     ,Rcpp::NumericVector par = Rcpp::NumericVector::create(2., 2., 1., 1., 0.)
                     ,Rcpp::NumericVector g_opt_par = Rcpp::NumericVector::create(.1, 1e-2, 100, 1e-3, 1e-5, 2)
                     ,Rcpp::NumericVector itv_opt_par = Rcpp::NumericVector::create(.01, 1e-3, 200, 1e-3, 1e-5, 5)
                     ){

  // check arguments
  int check_par_size = par.size();
  if(check_par_size != 5){
    Rcpp::stop("Number of default values for parameters must be equal to five.");
  }

  /* initial values */
  /* -------------- */

  // Name of parameters
  Rcpp::CharacterVector param_names =
    Rcpp::CharacterVector::create("bl", "br", "al", "ar", "m");

  /* store possibly provided values for parameters */
  unsigned is_m_provided = 0;

  // define optimization parameters
  struct multimin_params global_oparams =
    { (double)       g_opt_par[0]
     ,(double)       g_opt_par[1]
     ,(unsigned int) g_opt_par[2]
     ,(double)       g_opt_par[3]
     ,(double)       g_opt_par[4]
     ,(unsigned int) g_opt_par[5]
  };

  struct multimin_params interv_oparams =
    { (double)       itv_opt_par[0]
     ,(double)       itv_opt_par[1]
     ,(unsigned int) itv_opt_par[2]
     ,(double)       itv_opt_par[3]
     ,(double)       itv_opt_par[4]
     ,(unsigned int) itv_opt_par[5]
  };

  // casting the null value
  // necessary according to
  // https://gallery.rcpp.org/articles/optional-null-function-arguments/

  // initial value for m parameter
  if(provided_m_.isNotNull()){

    is_m_provided = 1;

    // casting the null value
    // necessary according to
    // https://gallery.rcpp.org/articles/optional-null-function-arguments/
    Rcpp::NumericVector provided_m(provided_m_);
    par[4] = provided_m[0];


  }

  /* store data */
  unsigned int Size = data.size();/*the number of data*/

  /* store guess */
  double   fmin=0;

  /* customized admin. of errors */
  /* --------------------------- */
  gsl_set_error_handler_off ();

  /* sort data */
  /* --------- */
  sortRcpp(data);

  /* initial values */
  /* -------------- */
  subboa_objf(data, 5, par, NULL, &fmin);

  /* output of initial values */
  /* ------------------------ */
  if(verb > 0){
    Rprintf("INITIAL VALUES\n");
    Rprintf("#  par    bl    br    al    ar    m     ll\n");
    Rprintf("#  value  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f\n", par[0], par[1], par[2], par[3], par[4],  fmin);
    Rprintf("\n");
  }

  /* no minimization */
  /* --------------- */
  if(method == 0){

    // main dataframe with coefficients
    Rcpp::DataFrame dt =
      Rcpp::DataFrame::create( Rcpp::Named("param") = param_names
                               ,Rcpp::Named("coef") = par
                               );

    Rcpp::List ans = Rcpp::List::create(
                                        Rcpp::Named("dt")              = dt
                                        ,Rcpp::Named("log-likelihood") = fmin
                                        );
    return ans;

  }


  /* global optimization  */
  /* -------------------- */
  if(method >= 1){

    // notice: par is updated by reference
    Rcpp::List g_opt_results =
      global_optim(
                   data
                   ,fmin
                   ,global_oparams
                   ,par
                   ,(unsigned) 5
                   ,(unsigned) 4
                   ,subboa_objf
                   ,subboa_objdf
                   ,subboa_objfdf
                   ,provided_m_
                   ,verb
                   );

    /* interval optimization  */
    /* ---------------------- */

    if(method >= 2){

      g_opt_results =
        interval_optim(
                       data
                       ,g_opt_results["type"]
                       ,g_opt_results["xmin"]
                       ,g_opt_results["xmax"]
                       ,par
                       ,g_opt_results["fmin"]
                       ,interv_oparams /* interval optimization parameters */
                       ,interv_step /* interval optimization step */
                       ,(unsigned) 5
                       ,(unsigned) 4
                       ,subboa_objf
                       ,subboa_objdf
                       ,subboa_objfdf
                       ,verb
                       );

      // updates log-likelihood
      fmin = Rcpp::as<double>(g_opt_results["fmin"]);
    }
  }

  // generate outputs

  // variables
  const size_t dim=(is_m_provided?4:5); /* set the size of the var-covar matrix */

  /* allocate var-covar matrix */
  // The V matrix has on its main diagonal the variance of parameters
  // on the lower diagonal the correlation coefficients
  // on the upper diagonal the covariances
  RcppGSL::Matrix V =  subboa_varcovar(par, Size, dim);
  // this matrix in its upper diagonal presents the covariances
  // and on its lower diagonal presents the correlation coefficients between the parameters

  // vector of standard errors
  Rcpp::NumericVector std_error =
    Rcpp::NumericVector::create(
                                 sqrt(V(0,0)) // bl
                                ,sqrt(V(1,1)) // br
                                ,sqrt(V(2,2)) // al
                                ,sqrt(V(3,3)) // ar
                                ,sqrt(V(4,4)) // m
                                );

  // main dataframe with coefficients
  Rcpp::List dt =
    Rcpp::DataFrame::create( Rcpp::Named("param")     = param_names
                            ,Rcpp::Named("coef")      = par
                            ,Rcpp::Named("std_error") = std_error
                      );

  // convert matrix
  Rcpp::NumericMatrix matrix = Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(V));

  // add names for the matrices
  colnames(matrix) = param_names;

  Rcpp::List ans = Rcpp::List::create(
                                       Rcpp::Named("dt")             = dt
                                      ,Rcpp::Named("log-likelihood") = fmin
                                      ,Rcpp::Named("matrix")         = matrix
                                     );


  return ans;

}
