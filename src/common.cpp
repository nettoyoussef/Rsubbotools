/*
  common.c (ver. 1.0.0) -- Various utilities functions for subbotools
  Copyright (C) 2003-2010 Giulio Bottazzi

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




double median(Rcpp::NumericVector data, size_t size){

  double median_value;

  // size is even
  if(size % 2 == 0){
    median_value = (data[size/2] + data[(size+1)/2])/2;
  }
  // size is odd
  else{
    median_value = data[(size+1)/2];
  }
  return median_value;
}

double calculate_index(size_t size){

  double index;

  // size is even
  if(size % 2 == 0){
    index = size/2;
  }
  // size is odd
  else{
    index = (size+1)/2;
  }
  return index;
}





void *my_alloc(size_t size){
  void *temp;
  if(!(temp = malloc(size))){
    Rcpp::stop("malloc, memory allocation failed");
  };
  return(temp);
}

void *my_realloc(void *ptr,size_t size){
  void *temp=NULL;
  if(size == 0){
    free(ptr);
  }
  else if(!(temp = realloc(ptr,size))){
    Rcpp::stop("realloc, memory allocation failed");
  };
  return(temp);
}

int sort_by_value(const void *a, const void *b){
  const double *da = (const double *) a;
  const double *db = (const double *) b;

  return (*da > *db) - (*da < *db);
}


// Conversion functions
// inspired by
// https://github.com/cran/Rmixmod/blob/master/src/Conversion.cpp
double * RcppNum_to_double(Rcpp::NumericVector x){

  int size = x.size();
  double *z = new double[size];

  for(int i =0; i<size; ++i){

    z[i] = x[i];

  }

  return z;

}

// Utilities

// [[Rcpp::export]]
void sortRcpp(Rcpp::NumericVector x){

  std::sort(x.begin(), x.end());

}


// Functions for optimization in the Subbotin Family

void set_m(
           Rcpp::NumericVector data
           ,Rcpp::NumericVector param
           ,unsigned m_position
           ,Rcpp::NumericVector xmin
           ,Rcpp::NumericVector xmax
           ,unsigned index
           ){

  // set boundaries on m
  xmin[m_position] = data[index];
  xmax[m_position] = data[index+1];

  // set the m value to the mid of the interval
  param[m_position] = .5*(xmin[m_position]+xmax[m_position]);

}



void check_new_minimum_core(
                            double *fmin
                            ,double dtmp1
                            ,Rcpp::NumericVector par
                            ,Rcpp::NumericVector xtmp
                            ){

  if(dtmp1 < *fmin){/* found new minimum */
    Rcppdeepcopy(xtmp, par); // set values to the last best value
    *fmin  = dtmp1;
  }
}


// this updates the objective function value,
// the parameters and the index position which
// references the m minimum
void check_new_minimum(
                       Rcpp::NumericVector data
                       ,double *fmin
                       ,double dtmp1
                       ,Rcpp::NumericVector par
                       ,Rcpp::NumericVector xtmp
                       ,unsigned i
                       ,unsigned *index
                       ,int verb
                       ){

  if(dtmp1 < *fmin){/* found new minimum */
    Rcppdeepcopy(xtmp, par); // set values to the last best value
    *fmin  = dtmp1;
    *index = i;

    // print results - we use arrows to
    // differentiate the global minimum
    if(verb > 1){

      int n = par.size();

      Rprintf("#>>> [%+.3e:%+.3e]", data[i], data[i+1]);
      for(int tmp=0; tmp < n;tmp++){
        Rprintf(" par[%d]= %e", tmp, par[i]);
      }
      Rprintf(" ll= %e\n\n", *fmin);
    }
  }else{
    // print results
    if(verb > 1){
      Rprintf("#    [%+.3e:%+.3e] ll=%e\n",
              data[i], data[i+1], dtmp1);
    }
  }
}

void Rcppdeepcopy(Rcpp::NumericVector x_orig, Rcpp::NumericVector x_dest){


  unsigned xo_size, xd_size;

  xo_size = x_orig.size();
  xd_size = x_dest.size();

  if(xo_size != xd_size){
    Rcpp::stop("Vectors must have the same size!");
  }

  for(int i =0; i<xo_size; ++i){
    x_dest[i] = x_orig[i];
  }
}

void create_parameter_names(char *par_names[], int size){

  if(size == 2){
    par_names[0] = "b";
    par_names[1] = "m";
  }
  else if(size == 4){
    par_names[0] = "bl";
    par_names[1] = "br";
    par_names[2] = "a";
    par_names[3] = "m";

  }
  else if(size == 5){
    par_names[0] = "bl";
    par_names[1] = "br";
    par_names[2] = "al";
    par_names[3] = "ar";
    par_names[4] = "m";
  }
}

void print_parameters(
                      Rcpp::NumericVector xtmp
                      ,Rcpp::IntegerVector type
                      ,Rcpp::NumericVector xmin
                      ,Rcpp::NumericVector xmax
                      ){

  int size = xtmp.size();
  char *par_names[size];

  create_parameter_names(par_names, size);

  Rprintf("#  par  value  type   xmin   xmax\n");
  for(int i=0; i<size; ++i){
    Rprintf("#  %s    %.2f      %i   %.2f   %.2f\n", par_names[i], xtmp[i], type[i], xmin[i], xmax[i]);
  }
  Rprintf("\n");
}

void print_results(
                   Rcpp::NumericVector par
                   ,double fmin
                   ){

  int size = par.size();
  char *par_names[size];

  create_parameter_names(par_names, size);

  Rprintf("#  par    ");
  for(int i =0; i<size; ++i){
    Rprintf("%s      ", par_names[i]);
  }
  Rprintf("ll\n");

  Rprintf("#  value  ");
  for(int i =0; i<size; ++i){
    Rprintf("%.3f  ", par[i]);
  }
  Rprintf("%.3f\n", fmin);
  Rprintf("\n");
}


// global optimization
Rcpp::List global_optim(
                        Rcpp::NumericVector data
                        ,double fmin
                        ,struct multimin_params global_oparams
                        ,Rcpp::NumericVector par
                        ,unsigned n_param
                        ,unsigned m_position
                        ,void (*f)    (Rcpp::NumericVector
                                       , const size_t
                                       , Rcpp::NumericVector, void *, double *)
                        ,void (* df)  (Rcpp::NumericVector
                                       , const size_t
                                       , Rcpp::NumericVector
                                       , void *
                                       , Rcpp::NumericVector)
                        ,void (* fdf) (Rcpp::NumericVector
                                       , const size_t
                                       , Rcpp::NumericVector
                                       , void *
                                       , double *
                                       , Rcpp::NumericVector)
                        ,Rcpp::Nullable<Rcpp::NumericVector> provided_m_ = R_NilValue
                        ,int verb = 0
                        ){


  /* ML estimation: global maximization */
  /* ---------------------------------- */

  // this function checks if the user requested the interruption of the
  // computation
  Rcpp::checkUserInterrupt();

  // declare variables
  Rcpp::IntegerVector type(n_param);
  Rcpp::NumericVector xmin(n_param);
  Rcpp::NumericVector xmax(n_param);

  // the initial guess take the value of the initial
  // parameters
  Rcpp::NumericVector xtmp(n_param); // xtmp[0] = b xtmp[1] = m
  Rcppdeepcopy(par, xtmp); // receive previous optimization value as starting point

  // creates temp log-likelihood
  double dtmp1;

  if(verb > 1){
    Rprintf("START OF GLOBAL OPTIMIZATION\n");
  }

  /* set initial minimization boundaries */
  /* ----------------------------------- */

  // set all parameters to be bounded in the interval
  // (a, +inf)
  type.fill(4);  // interval: (a,+inf); transf: x=a + \exp(y)
  xmin.fill(0);  // b must be positive
  xmax.fill(0);  // not used


  // correct the optimization criteria for m
  if(provided_m_.isNotNull()){

    // casting the null value
    // necessary according to
    // https://gallery.rcpp.org/articles/optional-null-function-arguments/
    Rcpp::NumericVector provided_m(provided_m_);

    // parameters for m
    type[m_position] = 3; //interval: [a,b]; transf: x= \frac{ a(1-\sin(y)) + b(1+\sin(y)}{2}
    xmin[m_position] = provided_m[0]; // uses the provided m to fix a value
    xmax[m_position] = provided_m[0]; // uses the provided m to fix a value
    xtmp[m_position] = provided_m[0];

    // for subbolafit, the transformation for b includes zero
    // when m is provided. Why?
    if(n_param == 4){

      // parameters for bl
      type[0] = 1;

      // parameters for br
      type[1] = 1;

    }
  }
  else{

    // parameters for m
    type[m_position] = 0; // interval: (-inf,+inf); transf: x=y
    xmin[m_position] = 0; // not used
    xmax[m_position] = 0; // not used
    xtmp[m_position] = par[m_position]; // receives value from previous step

  }

  if(verb > 1){
    Rprintf("Parameters for global optimization: \n");
    print_parameters(xtmp, type, xmin, xmax);
  }


  /* perform global minimization */
  /* --------------------------- */

  multimin(data                   // sample
           ,n_param               // number of parameters
           ,xtmp                  // starting guess / returns parameters value
           ,&dtmp1                // pointer to update the minimum likelihood
           ,type                  // type of transformation of the data
           ,xmin                  // minimum values for the parameters
           ,xmax                  // maximum values for the parameters
           ,f                     // objective function to minimize
           ,df                    // df/dx of the objective function
           ,fdf                   // objf + objdf
           ,NULL                  // fparams
           ,global_oparams        // parameters for the optimization
           ,verb                  // set verbosity level
           );

  // updates minimum values and parameters
  check_new_minimum_core(&fmin, dtmp1, par, xtmp);

  if(verb > 1){
    Rprintf("Results of global optimization: \n");
    print_results(par, fmin);
    Rprintf("END OF GLOBAL OPTIMIZATION\n");
  }

  Rcpp::List ans =
    Rcpp::List::create(
                       // this has the same address as the parent object
                       // I only bring it here because I have to bring the other values
                       // as well
                       // Rcpp::Named("par")  = par
                       Rcpp::Named("type") = Rcpp::wrap(type)
                       ,Rcpp::Named("xmin") = Rcpp::wrap(xmin)
                       ,Rcpp::Named("xmax") = Rcpp::wrap(xmax)
                       ,Rcpp::Named("fmin") = fmin
                       );

  return ans;
}




// - interv_step - interval optimization parameters
// - interv_oparams - interval optimization parameters
Rcpp::List interval_optim(
                          Rcpp::NumericVector data
                          ,Rcpp::IntegerVector type
                          ,Rcpp::NumericVector xmin
                          ,Rcpp::NumericVector xmax
                          ,Rcpp::NumericVector par
                          ,double fmin
                          ,struct multimin_params interv_oparams
                          ,int interv_step
                          ,unsigned n_param
                          ,unsigned m_position
                          ,void (*f)    (Rcpp::NumericVector
                                         , const size_t
                                         , Rcpp::NumericVector
                                         , void *
                                         , double *)
                          ,void (* df)  (Rcpp::NumericVector
                                         , const size_t
                                         , Rcpp::NumericVector
                                         , void *
                                         , Rcpp::NumericVector)
                          ,void (* fdf) (Rcpp::NumericVector
                                         , const size_t
                                         , Rcpp::NumericVector
                                         , void *
                                         , double *
                                         , Rcpp::NumericVector)
                          ,int verb
                          ){

  // local
  double dtmp1;
  Rcpp::NumericVector xtmp(n_param); /* store the temporary minimum  */
  unsigned int index;
  unsigned int oldindex;
  unsigned int i;
  unsigned int max;
  unsigned int min;
  unsigned int size = data.size();
  //unsigned int n = 2; // number of parameters

  // only x[0] = b = par[0] and x[1]=m = par[2] are optimized

  // for Subbofit, you pass only b and m
  // for subboafit, you pass all parameters
  // for subbolafit, you pass all parameters (?)
  // necessary to pass value for the b parameter
  Rcppdeepcopy(par, xtmp);

  /* perform interval minimization */
  /* ----------------------------- */

  if(verb > 1){
    Rprintf("INTERVAL-WISE OPTIMIZATION\n");
  }


  /* check plausibility of specified number of intervals */
  if(interv_step>size/2){
    Rprintf("#WARNING: Too much intervals specified. Number of intervals set to %d\n"
            ,size/2);
    interv_step = size/2;
  }

  // m is always bounded by xmin and xmax
  type[m_position] = 3;

  // this finds the index of the maximum value below
  // the estimated location parameter
  // the loop advances until it finds the first data
  // point that is bigger than m
  for(i = 0; data[i]<=par[m_position] && i<size; i++);

  // this corrects the value of the index
  // to get the value imediatelly below m
  if(i == 0) // this corrects for 2 data points
    index = 0;
  // not sure why this exists, can't think of a use case
  else if (i == size)
    index= size-2;
  // default route, it becomes one data point below
  // the m value
  else index = i-1;

  // set m value
  set_m(data, xtmp, m_position, xmin, xmax, index);

  if(verb > 1){
    Rprintf("INITIAL INTERVAL-WISE OPTIMIZATION\n");
  }

  // initial optimization
  // minimization
  multimin(
           data             // sample
           ,n_param         // number of parameters
           ,xtmp            // starting guess / returns parameters value
           ,&dtmp1          // pointer to update the minimum likelihood
           ,type            // type of transformation of the data
           ,xmin            // minimum values for the parameters
           ,xmax            // maximum values for the parameters
           ,f               // objective function to minimize
           ,df              // df/dx of the objective function
           ,fdf             // objf + objdf
           ,NULL            // fparams
           ,interv_oparams  // parameters for the optmization
           ,verb            // set verbosity level
           );

  // updates minimum values and parameters
  // we use arrows #>>> to represent a new minimum
  check_new_minimum(
                    data
                    ,&fmin
                    ,dtmp1
                    ,par
                    ,xtmp
                    ,index  // this is not used here, we just maintain the same index
                    ,&index
                    ,verb
                    );


  // set values to the last best value
  //Rcppdeepcopy(par, xtmp)
  // to check: it seems that the original code doesn't update the
  // 'ar' value for subboafit

  // this puts the index in the value that is just below the
  // estimate for the m parameter
  max = min = index;


  if(verb > 1){
    Rprintf("RIGHT-SIDE OPTIMIZATION\n");
  }

  /* move to the right compute new local minima and compare with
     global minimum */
  do{
    oldindex = index;

    // this function checks if the user requested the interruption of the
    // computation
    Rcpp::checkUserInterrupt();

    // the algorithm moves interv_step (default:10) steps to the right to see
    // if it finds a new minimum
    // a step is a position on the vector of data
    //  second condition i<size-1 avoids steeeping out of vector strip
    for(i = max+1; i<=max+interv_step && i<size-1; i++){

      // set initial condition
      Rcppdeepcopy(par, xtmp); // set values to the last best value

      // set m value
      set_m(data, xtmp, m_position, xmin, xmax, i);

      // initial optimization
      // minimization
      multimin(
               data             // sample
               ,n_param         // number of parameters
               ,xtmp            // starting guess / returns parameters value
               ,&dtmp1          // pointer to update the minimum likelihood
               ,type            // type of transformation of the data
               ,xmin            // minimum values for the parameters
               ,xmax            // maximum values for the parameters
               ,f               // objective function to minimize
               ,df              // df/dx of the objective function
               ,fdf             // objf + objdf
               ,NULL            // fparams
               ,interv_oparams  // parameters for the optmization
               ,verb            // set verbosity level
               );

      // updates minimum values and parameters
      check_new_minimum(
                        data
                        ,&fmin
                        ,dtmp1
                        ,par
                        ,xtmp
                        ,i
                        ,&index
                        ,verb
                        );

    }
    // counts the number of intervals explored
    // and adjusts max for the next iteration, if it occurs
    max = i-1;
  }
  while(index!=oldindex);

  if(verb > 1){
    Rprintf("LEFT-SIDE OPTIMIZATION\n");
  }

  // move to the left compute new local minima and compare with
  // best local minimum
  do{
    oldindex = index;

    // this function checks if the user requested the interruption of the
    // computation
    Rcpp::checkUserInterrupt();

    for(i = min-1; (int) i >= (int) min-interv_step && (int) i >= 0; i--){

      // set initial condition
      Rcppdeepcopy(par, xtmp); // set values to the last best value

      // set m value
      set_m(data, xtmp, m_position, xmin, xmax, i);

      // initial optimization
      // minimization
      multimin(
               data             // sample
               ,n_param         // number of parameters
               ,xtmp            // starting guess / returns parameters value
               ,&dtmp1          // pointer to update the minimum likelihood
               ,type            // type of transformation of the data
               ,xmin            // minimum values for the parameters
               ,xmax            // maximum values for the parameters
               ,f               // objective function to minimize
               ,df              // df/dx of the objective function
               ,fdf             // objf + objdf
               ,NULL            // fparams
               ,interv_oparams  // parameters for the optmization
               ,verb            // set verbosity level
               );

      // updates minimum values and parameters
      check_new_minimum(
                        data
                        ,&fmin
                        ,dtmp1
                        ,par
                        ,xtmp
                        ,i
                        ,&index
                        ,verb
                        );
    }
    min = i+1;
  }
  while(index!=oldindex);

  /* store final values */
  /* ------------------ */

  if(verb > 1){
    Rprintf("Results of interval optimization: \n");
    print_results(par, fmin);
    Rprintf("\n");
    Rprintf("#  intervals explored: %d\n",max-min);
    Rprintf("\n");
  }

  Rcpp::List ans =
    Rcpp::List::create(
                       Rcpp::Named("fmin") = fmin
                       );

  return ans;
}





/* functions for the information matrix of asymmetric subbotin  */
/* ----------------------------------------------------------- */

double B0(double x){

  return pow(x,1./x)*gsl_sf_gamma(1.+1./x);

}


double B1(double x){
  const double dt1 = 1.+1./x;
  return pow(x,1./x-1.)*gsl_sf_gamma(dt1)*(log(x)+gsl_sf_psi(dt1));

}


double B2(double x){
  const double dt1 = 1.+1./x;
  const double dt2 = log(x);
  const double dt3 = gsl_sf_psi(dt1);
  const double dt4 = gsl_sf_psi_1(dt1);

  return pow(x,1./x-2.)*gsl_sf_gamma(dt1)*(dt2*dt2+2*dt2*dt3+dt3*dt3+dt4);

}


double dB0dx(double x){

  return B0(x)/(x*x)-B1(x)/x;

}

double dB0dx2(double x){

  const double dt1=x*x;

  return -B0(x)*(3.*x-1.)/(dt1*dt1) + 2.*B1(x)*(x-1.)/(dt1*x) + B2(x)/dt1;

}



/* long options */
/* ------------ */

struct option gb_long_options[] = {
  {"version", no_argument,       NULL,  0 },
  {"help", no_argument,       NULL,  'h' },
  {0,         0,                 0,  0 }
};

int gb_option_index = 0;
