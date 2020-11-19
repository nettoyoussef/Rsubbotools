/*
  multimin.c (ver. 1.2) -- Interface to GSL multidim. minimization
  Copyright (C) 2002-2014 Giulio Bottazzi

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  (version 2) as published by the Free Software Foundation;
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program. If not see <http://www.gnu.org/licenses/>.

*/ 


/*

multimin is an interface to the various GSL minimization
routines. When invoked, all the information necessary to perform the
minimization are passed as formal parameters. This generate a pretty
long, FORTRAN-like, list of parameters. This approach allows, however,
to black-box, as far as possible, the interior functioning of the
routine from the rest of the program.

Let's analyse the calling convention in details:

multimin(size_t n,double *x,double *fun,
         unsigned *type,double *xmin,double *xmax,
         void (*f) (const size_t,const double *,void *,double *),
         void (* df) (const size_t,const double *, void *,double *),
         void (* fdf) (const size_t,const double *, void *,double *,double *),
         void *fparams,
         const struct multimin_params oparams)

where

--------------------------------------------------------------------
n

INPUT: dimension of the problem, number of independent variables of
the function.

--------------------------------------------------------------------
x

INPUT: pointer to an array of n values x[0],...x[n-1] containing the
initial estimate of the minimum point
OUTPUT: contains the final estimation of the minimum position

--------------------------------------------------------------------
type

a pointer to an array of integer type[1],...,type[n-1] describing the
boundary conditions for the different variables. The problem is solved
as an unconstrained one on a suitably transformed variable y. Possible
values are:

  Interval:                                       Transformation:
  0 unconstrained                                 x=y
  1 semi-closed right half line [ xmin,+infty )   x=xmin+y^2
  2 semi-closed left  half line ( -infty,xmax ]   x=xmax-y^2
  3 closed interval              [ xmin,xmax ]    x=SS+SD*sin(y)
  4 open right half line        ( xmin,+infty )   x=xmin+exp(y)
  5 open left  half line        ( -infty,xmax )   x=xmax-exp(y)
  6 open interval                ( xmin,xmax )    x=SS+SD*tanh(y)

where SS=.5(xmin+xmax) SD=.5(xmax-xmin)

There are also other UNSUPPORTED transformations used in various test
  7 open interval                ( xmin,xmax )    x=SS+SD*(1+y/sqrt(1+y^2))
  8 open right half line        ( xmin,+infty )   x=xmin+.5*(y+sqrt(1+y^2))
  9 open left  half line        ( -infty,xmax )   x=xmax+.5*(y-sqrt(1+y^2))
--------------------------------------------------------------------
xmin 
xmax

pointers to arrays of double containing respectively the lower and
upper boundaries of the different variables. For a given variable,
only the values that are implied by the type of constraints, defined
as in *type, are actually inspected.

--------------------------------------------------------------------
f                

f calculates the objective function at a specified point x. Its
specification is

void (*f) (const size_t n, const double *x,void *fparams,double *fval)

      n
      INPUT: the number of variables

      x
      INPUT:the point at which the function is required

      fparams
      pointer to a structure containing parameters required by the
      function. If no external parameter are required it can be set to
      NULL.

      fval 
      OUTPUT: the value of the objective function at the current point
      x.

--------------------------------------------------------------------
df               

df calculates the gradient of the objective function at a specified
point x. Its specification is

void (*df) (const size_t n, const double *x,void *fparams,double *grad)

      n
      INPUT: the number of variables

      x
      INPUT:the point at which the function is required

      fparams
      pointer to a structure containing parameters required by the
      function. If no external parameter are required it can be set to
      NULL.

      grad
      OUTPUT: the values of the gradient of the objective function at
      the current point x are stored in grad[0],...,grad[n-1].

--------------------------------------------------------------------
fdf              

fdf calculates the value and the gradient of the objective function at
a specified point x. Its specification is

void (*fdf) (const size_t n, const double *x,void *fparams,double *fval,double *grad)

      n
      INPUT: the number of variables

      x
      INPUT:the point at which the function is required

      fparams
      pointer to a structure containing parameters required by the
      function. If no external parameter are required it can be set to
      NULL.

      fval 
      OUTPUT: the value of the objective function at the current point
      x.

      grad
      OUTPUT: the values of the gradient of the objective function at
      the current point x are stored in grad[0],...,grad[n-1].

--------------------------------------------------------------------
fparams

pointer to a structure containing parameters required by the
function. If no external parameter are required it can be set to NULL.

--------------------------------------------------------------------

oparams

structure of the type "multimin_params" containing the optimization
parameters. The members are

      double step_size
          size of the first trial step 

      double tol
          accuracy of the line minimization

      unsigned maxiter
          maximum number of iterations 

      double epsabs
          accuracy of the minimization

      double maxsize;
          final size of the simplex

      unsigned method
          method to use. Possible values are:

          0: Fletcher-Reeves conjugate gradient
          1: Polak-Ribiere conjugate gradient
          2: Vector Broyden-Fletcher-Goldfarb-Shanno method
          3: Steepest descent algorithm
          4: Nelder-Mead simplex
          5: Vector Broyden-Fletcher-Goldfarb-Shanno method ver. 2
          6: Simplex algorithm of Nelder and Mead ver. 2
          7: Simplex algorithm of Nelder and Mead: random initialization

      unsigned verbosity
          if greater then 0 print info on intermediate steps

*/
#include "multimin.h"

struct g_params {

  // sample to estimate the parameters/variables (rows)
  Rcpp::NumericVector data;
  
  // dimension of the problem (number of parameters/variables) (columns)
  size_t n;  
  
  // which transformation to do for each parameter/variable, if any
  Rcpp::IntegerVector type;

  // vector with minimum value for each parameter
  Rcpp::NumericVector xmin;

  // vector with maximum value for each parameter
  Rcpp::NumericVector xmax;

  // depending on the algorithm, one of the following
  // functions is used for the minimization process
  
  // *f - function to be minimized
  // Arguments, on this order, are
  // vector of data
  // size of vector
  // initial guess pointer - is substituted by the final solution's position
  // parameters to the minimization problem
  // solution - minimum value obtained
  void (*f) (Rcpp::NumericVector, const size_t, std::vector<double>, void *, double *);
    
  // *df - derivative of the function to be minimized
  // Arguments, on this order, are
  // vector of data
  // size of vector
  // initial guess pointer - is substituted by the final solution's position
  // parameters to the minimization problem
  // solution - minimum value obtained
  void (* df) (Rcpp::NumericVector, const size_t, std::vector<double>, void *, std::vector<double>);
    
  // combination of f and df
  void (* fdf) (Rcpp::NumericVector, const size_t, std::vector<double>, void *, double *, std::vector<double>);

  // parameters for the functions
  void *fparams;
};


void do_data_transformation( std::vector<double> x
                            ,const size_t n
                            ,const gsl_vector *y
                            ,Rcpp::IntegerVector type
                            ,Rcpp::NumericVector xmin
                            ,Rcpp::NumericVector xmax
                            ){

  
  size_t i;
  double dtmp1;

  for(i=0;i<n;i++){
    if(type[i] == NA_INTEGER){
      x[i]= gsl_vector_get(y,i);
    }
    else
      switch(type[i]){
      case 0:/* (-inf,+inf) */
        x[i]= gsl_vector_get(y,i); 
        break;
      case 1:/* [a,+inf) */
        x[i]= xmin[i]+gsl_vector_get(y,i)*gsl_vector_get(y,i); 
        break;
      case 2:/* (-inf,a] */
        x[i]= xmax[i]-gsl_vector_get(y,i)*gsl_vector_get(y,i); 
        break;
      case 3:/* [a,b] */
        dtmp1 = sin( gsl_vector_get(y,i) );
        x[i]= .5*(xmin[i]*(1-dtmp1) +xmax[i]*(1+dtmp1)); 
        break;
      case 4:/* (a,+inf) */
        dtmp1 = exp( gsl_vector_get(y,i) );
        x[i]= xmin[i]+dtmp1;
        break;
      case 5:/* (-inf,a) */
        dtmp1 = -exp( gsl_vector_get(y,i) );
        x[i]= xmax[i]+dtmp1;
        break;
      case 6:/* (a,b) */
        dtmp1 = tanh( gsl_vector_get(y,i) );
        x[i]= .5*(xmin[i]*(1-dtmp1) +xmax[i]*(1+dtmp1));
        break;
      case 7:/* (a,b) second approach */
        dtmp1 = gsl_vector_get(y,i)/sqrt(1.+gsl_vector_get(y,i)*gsl_vector_get(y,i));
        x[i]= .5*(xmin[i]*(1-dtmp1) +xmax[i]*(1+dtmp1));
        break;
      case 8:/* (a,+inf) second approach */
        dtmp1 = sqrt(1.+gsl_vector_get(y,i)*gsl_vector_get(y,i));
        x[i]= xmin[i] + .5*(gsl_vector_get(y,i)+dtmp1);
        break;
      case 9:/* (-inf,a) second approach */
        dtmp1 = sqrt(1.+gsl_vector_get(y,i)*gsl_vector_get(y,i));
        x[i]= xmax[i] + .5*(gsl_vector_get(y,i)-dtmp1);
        break;
      }
  }
}
  
  
void do_data_transformation_df( std::vector<double> x
                               ,const size_t n
                               ,const gsl_vector *y
                               ,Rcpp::IntegerVector type
                               ,Rcpp::NumericVector xmin
                               ,Rcpp::NumericVector xmax
                               ,std::vector<double> dx
                               ){

  size_t i;
  double dtmp1,dtmp2;
  
  for(i=0;i<n;i++){
    if(type[i] == NA_INTEGER){
      dx[i]= 1;
    }
    else
      switch(type[i]){
      case 0:/* (-inf,+inf) */
        dx[i]= 1;
        break;
      case 1:/* [a,+inf) */
        dx[i]= 2.*gsl_vector_get(y,i);
        break;
      case 2:/* (-inf,a] */
        dx[i]= -2.*gsl_vector_get(y,i);
        break;
      case 3:/* [a,b] */
        dx[i]= .5*(xmax[i]-xmin[i])*cos(gsl_vector_get(y,i));
        break;
      case 4:/* (a,+inf) */
        dx[i]= exp( gsl_vector_get(y,i) );
        break;
      case 5:/* (-inf,a) */
        dx[i]= -exp( gsl_vector_get(y,i) );
        break;
      case 6:/* (a,b) */
        dtmp1 = cosh( gsl_vector_get(y,i) );
        dx[i]= .5*(xmax[i]-xmin[i])/(dtmp1*dtmp1);
        break;
      case 7:/* (a,b) second approach */
        dtmp1 = (1.+gsl_vector_get(y,i)*gsl_vector_get(y,i))*sqrt(1.+gsl_vector_get(y,i)*gsl_vector_get(y,i)) ;
        dx[i]= .5*(xmax[i]-xmin[i])/dtmp1;
        break;
      case 8:/* (a,+inf) second approach */
        dtmp1 = gsl_vector_get(y,i);
        dtmp2 = sqrt(1.+dtmp1*dtmp1);
        dx[i]= .5*(dtmp1+dtmp2)/dtmp2;
        break;
      case 9:/* (-inf,a) second approach */
        dtmp1 = gsl_vector_get(y,i);
        dtmp2 = sqrt(1.+dtmp1*dtmp1);
        dx[i]= .5*(dtmp2-dtmp1)/dtmp2;
        break;
      }
  }
 } 


// this function returns the evaluation of f(x_1, x_2, x_n)
// for a guess y_1, ..., y_n
// in this context, y represents the coefficients to be
// estimated
static double g(const gsl_vector *y, void *gparams){

  struct g_params *p= (struct g_params *) gparams;
  
  double res=GSL_NAN;/* the function is forced to return a value */
  

  /* dereference useful stuff */
  const size_t n = p->n;
  Rcpp::IntegerVector type = p->type;
  Rcpp::NumericVector xmin = p->xmin;
  Rcpp::NumericVector xmax = p->xmax;

  //double *x = (double *) multimin_alloc(sizeof(double)*n);
  std::vector<double> x(n);

  // variable transformation
  
  // type is used to select the boundary conditions for
  // the values of the different variables/coefficents.
  // the process uses a data transformation to convert a
  // constrained problem in a unconstrained one.

  // Notice: type is also a vector, so for each 
  // variable a different transformation can be selected
  do_data_transformation(x, n, y, type, xmin, xmax);

  p->f(p->data, n,x,p->fparams,&res);
  //free (x);
  return(res);

}



// this function must allocate in a pointer the
// value of the gradient of f evaluated at point
// y_1, ..., y_n
static void dg(const gsl_vector *y, void *gparams, gsl_vector *dg){

  struct g_params *p= (struct g_params *) gparams;
  
  
  /* dereference useful stuff */
  const size_t n = p->n;
  Rcpp::IntegerVector type = p->type;
  Rcpp::NumericVector xmin = p->xmin;
  Rcpp::NumericVector xmax = p->xmax;

  // x is the guess after the variable transformation
  std::vector<double> x(n);
  //double *x = (double *) multimin_alloc(sizeof(double)*n);

  // dx is the gradient after the variable transformation
  std::vector<double> dx(n);
  //double *dx = (double *) multimin_alloc(sizeof(double)*n);

  std::vector<double> df(n);
  //double *df = (double *) multimin_alloc(sizeof(double)*n);
  

  // variable transformation
  
  // type is used to select the boundary conditions for
  // the values of the different variables/coefficents.
  // the process uses a data transformation to convert a
  // constrained problem in a unconstrained one.

  // Notice: type is also a vector, so for each 
  // variable a different transformation can be selected
  
  // data transf. for x  
  do_data_transformation(x, n, y, type, xmin, xmax);

  // data transf. for dx
  do_data_transformation_df(x, n, y, type, xmin, xmax, dx);

  // find the gradient at x and update df with it 
  p->df(p->data, n, x, p->fparams, df);

  /* debug output; comment out if necessary */
  /*   fprintf(stderr,"#dg: x=( "); */
  /*   for(i=0;i<n;i++) */
  /*     fprintf(stderr,"%f ",gsl_vector_get(x,i)); */
  /*   fprintf(stderr,") dx=( "); */
  /*   for(i=0;i<n;i++) */
  /*     fprintf(stderr,"%f ",gsl_vector_get(dx,i)); */
  /*   fprintf(stderr,") df=( "); */
  /*   for(i=0;i<n;i++) */
  /*     fprintf(stderr,"%f ",gsl_vector_get(dg,i)); */

  // apply the data transformation to gradient df with the dx vector  
  // and set it to dg
  for(unsigned i=0; i<n;i++){
    gsl_vector_set(dg,i,df[i]*dx[i]);
  }

  /* debug output; comment out if necessary */
  /*   fprintf(stderr,") dg=( "); */
  /*   for(i=0;i<n;i++) */
  /*     fprintf(stderr,"%f ",gsl_vector_get(dg,i)); */
  /*   fprintf(stderr,")\n"); */

  //free (x);
  //free (dx);
  //free (df);

}

// this function combines g and dg, and the gsl package chooses
// which to evaluate or both, depending on the algorithm used
static void gdg(const gsl_vector *y, void *gparams, double *g, gsl_vector *dg){

  struct g_params *p= (struct g_params *) gparams;
  
  size_t i;

  /* dereference useful stuff */
  const size_t n = p->n;
  Rcpp::IntegerVector type = p->type;
  Rcpp::NumericVector xmin = p->xmin;
  Rcpp::NumericVector xmax = p->xmax;

  std::vector<double> x(n);
  std::vector<double> dx(n);
  std::vector<double> df(n);
  
  //double *x  = (double *) multimin_alloc(sizeof(double)*n);
  //double *dx = (double *) multimin_alloc(sizeof(double)*n);
  //double *df = (double *) multimin_alloc(sizeof(double)*n);

  // variable transformation
  
  // type is used to select the boundary conditions for
  // the values of the different variables/coefficents.
  // the process uses a data transformation to convert a
  // constrained problem in a unconstrained one.

  // Notice: type is also a vector, so for each 
  // variable a different transformation can be selected

  // data transf. for x  
  do_data_transformation(x, n, y, type, xmin, xmax);

  // data transf. for dx
  do_data_transformation_df(x, n, y, type, xmin, xmax, dx);

  // find the gradient at x and update df with it 
  p->fdf(p->data, n, x, p->fparams, g, df);

  /* debug output; comment out if necessary */    
  /*   fprintf(stderr,"#gdg: f=%f x=( ",g); */
  /*   for(i=0;i<n;i++) */
  /*     fprintf(stderr,"%f ",gsl_vector_get(x,i)); */
  /*   fprintf(stderr,") dx=( "); */
  /*   for(i=0;i<n;i++) */
  /*     fprintf(stderr,"%f ",gsl_vector_get(dx,i)); */
  /*   fprintf(stderr,") df=( "); */
  /*   for(i=0;i<n;i++) */
  /*     fprintf(stderr,"%f ",gsl_vector_get(dg,i)); */
  /*   fprintf(stderr,")\n"); */
  
  // apply the data transformation to gradient df with the dx vector  
  // and set it to dg
  for(i=0;i<n;i++){
    gsl_vector_set(dg,i,df[i]*dx[i]);
  }
  
  //free (x);
  //free (dx);
  //free (df);
}

// Select the algorithm to be used in the minimization process
// this function has to be this bizarre turnaround because
// gsl only accepts const arguments for its minimizers.
struct multimin_algorithm choose_algorithm(unsigned int method){


  Rprintf("Choosing algorithm:\n");
  
  const gsl_multimin_fdfminimizer_type *Tfdf;
  const gsl_multimin_fminimizer_type *Tf;
  const char *Tname;

  
  /* set the algorithm */
  switch(method){
  case 0:/* Fletcher-Reeves conjugate gradient */
    Tfdf = gsl_multimin_fdfminimizer_conjugate_fr;
    Tname = Tfdf->name;
    break;
  case 1:/* Polak-Ribiere conjugate gradient */
    Tfdf = gsl_multimin_fdfminimizer_conjugate_pr;
    Tname = Tfdf->name;
    break;
  case 2:/* Vector Broyden-Fletcher-Goldfarb-Shanno method */
    Tfdf = gsl_multimin_fdfminimizer_vector_bfgs;
    Tname = Tfdf->name;
    break;
  case 3:/* Steepest descent algorithm */
    Tfdf =gsl_multimin_fdfminimizer_steepest_descent;
    Tname = Tfdf->name;
    break;
  case 4:/* Simplex algorithm of Nelder and Mead */
    Tf = gsl_multimin_fminimizer_nmsimplex;
    Tname = Tf->name;
    break;
  case 5:/*  Vector Broyden-Fletcher-Goldfarb-Shanno2 method */
    Tfdf = gsl_multimin_fdfminimizer_vector_bfgs2;
    Tname = Tfdf->name;
    break;
  case 6:/* Simplex algorithm of Nelder and Mead version 2 */
    Tf = gsl_multimin_fminimizer_nmsimplex2;
    Tname = Tf->name;
    break;
  case 7:/* Simplex algorithm of Nelder and Mead: random initialization */
    Tf = gsl_multimin_fminimizer_nmsimplex2rand;
    Tname = Tf->name;
   break;
   
  default:
    Rprintf("Optimization method not recognized. Specify one of the following:\n\n");

    Rprintf("0: Fletcher-Reeves conjugate gradient\n");
    Rprintf("1: Polak-Ribiere conjugate gradient\n");
    Rprintf("2: Vector Broyden-Fletcher-Goldfarb-Shanno method\n");
    Rprintf("3: Steepest descent algorithm\n");
    Rprintf("4: Nelder-Mead simplex\n");
    Rprintf("5: Vector Broyden-Fletcher-Goldfarb-Shanno method ver. 2\n");
    Rprintf("6: Simplex algorithm of Nelder and Mead ver. 2\n");
    Rprintf("7: Simplex algorithm of Nelder and Mead: random initialization\n");

    Rcpp::stop("Choose again with one of the methods above.");
  }


  struct multimin_algorithm multimin_alg;
  
  multimin_alg.Tfdf  = Tfdf;
  multimin_alg.Tf    = Tf;
  multimin_alg.Tname = Tname;

  Rprintf("Algorithm chosen: %s\n", multimin_alg.Tname);

  return multimin_alg;
  
}


// compute values of y for initial condition
void do_initial_data_transform( std::vector<double> x
                               ,const size_t n
                               ,gsl_vector *y
                               ,Rcpp::IntegerVector type
                               ,Rcpp::NumericVector xmin
                               ,Rcpp::NumericVector xmax
                               ){


  size_t i;
  double dtmp1;


  Rprintf("DATA TRANSFORMATION\n");

  
  Rprintf("#    - variables initial value and boundaries\n");

  for(i=0;i<n;i++){
    if(type[i] == NA_INTEGER){
      gsl_vector_set(y,i,x[i]);
      Rprintf("#    x[%d]=%.2f (-inf,+inf)  trans 0 -> %.2f\n",(int) i,x[i],gsl_vector_get(y,i));
    }
    else
      switch(type[i]){
      case 0:/* (-inf,+inf) */
        gsl_vector_set(y,i,x[i]);
        Rprintf("#    x[%d]=%.2f (-inf,+inf)  trans 0 -> %.2f\n",(int) i,x[i],gsl_vector_get(y,i));
        break;
      case 1:/* [a,+inf) */
        gsl_vector_set(y,i,sqrt( x[i]-xmin[i] ));
        Rprintf("#    x[%d]=%.2f [%.3g,+inf)  trans 1 -> %.2f\n",(int) i,x[i],xmin[i],gsl_vector_get(y,i));
        break;
      case 2:/* (-inf,a] */
        gsl_vector_set(y,i,sqrt( xmax[i]-x[i] ));
        Rprintf("#    x[%d]=%.2f (-inf,%.3f]        trans 2 -> %.2f\n",(int) i,x[i],xmax[i],gsl_vector_get(y,i));
        break;
      case 3:/* [a,b] */
        dtmp1 = (xmax[i]>xmin[i]?
                 (2.*x[i]-xmax[i]-xmin[i])/(xmax[i]-xmin[i]) : 0);
        /*       dtmp1 = (2.*x[i]-xmax[i]-xmin[i])/(xmax[i]-xmin[i]); */
        gsl_vector_set(y,i,asin( dtmp1 ));
        Rprintf("#    x[%d]=%.2f [%.3f,%.3f] trans 3 -> %.2f\n",(int) i,x[i],xmin[i],xmax[i],gsl_vector_get(y,i));
        break;
      case 4:/* (a,+inf) */
        gsl_vector_set(y,i,log( x[i]-xmin[i] ));
        Rprintf("#    x[%d]=%.2f (%.3f,+inf)  trans 4 -> %.2f\n",(int) i,x[i],xmin[i],gsl_vector_get(y,i));
        break;
      case 5:/* (-inf,a) */
        gsl_vector_set(y,i,log( xmax[i]-x[i] ));
        Rprintf("#    x[%d]=%.2f (-inf,%.3f)  trans 5 -> %.2f\n",(int) i,x[i],xmax[i],gsl_vector_get(y,i));
        break;
      case 6:/* (a,b) */
        dtmp1 = (2.*x[i]-xmax[i]-xmin[i])/(xmax[i]-xmin[i]);
        gsl_vector_set(y,i,gsl_atanh ( dtmp1 ));
        Rprintf("#    x[%d]=%.2f (%.3f,%.3f) trans 6 -> %.2f\n",(int) i,x[i],xmin[i],xmax[i],gsl_vector_get(y,i));
        break;
      case 7:/* (a,b) second approach */
        dtmp1 = (2.*x[i]-xmax[i]-xmin[i])/(xmax[i]-xmin[i]);
        gsl_vector_set(y,i, dtmp1/sqrt(1-dtmp1*dtmp1));
        Rprintf("#    x[%d]=%.2f (%.3f,%.3f) trans 7 -> %.2f\n",(int) i,x[i],xmin[i],xmax[i],gsl_vector_get(y,i));
        break;
      case 8:/* (a,+inf) second approach */
        dtmp1 = x[i]-xmin[i];
        gsl_vector_set(y,i, dtmp1-1./(4.*dtmp1));
        Rprintf("#    x[%d]=%.2f (%.3f,+inf)  trans 8 -> %.2f\n",(int) i,x[i],xmin[i],gsl_vector_get(y,i));
        break;
      case 9:/* (-inf,a) second approach */
        dtmp1 = xmax[i]-x[i];
        gsl_vector_set(y,i, 1./(4.*dtmp1)-dtmp1);
        Rprintf("#    x[%d]=%.2f (-inf,%.3f)  trans 9 -> %.2f\n",(int) i,x[i],xmax[i],gsl_vector_get(y,i));
        break;
      }
  }


  Rprintf("END OF DATA TRANSFORMATION\n");
  
}




//' - n - INPUT: dimension of the problem, number of independent variables of the function. 
//' 
//' - x - INPUT: pointer to an array of n values x[0],...x[n-1] containing the initial
//'       estimate of the minimum point. OUTPUT: contains the final estimation of the
//'       minimum position 
//' 
//' - fmin - OUPUT: return the value of the function at the minimum 
//' 
//' - type - INPUT a pointer to an array of n integers type[1],...,type[n-1] describing the
//'          boundary conditions for the different variables. The problem is solved as an
//'          unconstrained one on a suitably transformed variable y. Possible values are:
//'          value interval      data transformation          obs
//'           0:   (-inf,+inf)   x=y
//'           1:   [a,+inf)      x=a+y^2
//'           2:   (-inf,a]      x=a-y^2
//'           3:   [a,b]         x= \frac{ a(1-\sin(y)) + b(1+\sin(y)}{2}
//'           4:   (a,+inf)      x=a + \exp(y)
//'           5:   (-inf,a)      x=a - \exp(y)
//'           6:   (a,b)         x= \frac{ a(1-\tanh(y)) + b(1+\tanh(y)}{2}
//'           7:   (a,b)         x=\frac{ a(1 - \frac{y}{ (1+y^2)^{1/2}} ) + b( 1 + \frac{y}{ (1+y^2)^{1/2}} )}{2}    second approach 
//'           8:   (a,+inf)      x= a + \frac{y + (1 + y^2)^{1/2} }{2}   second approach 
//'           9:   (-inf,a)      x=a + \frac{y - (1 + y^2)^{1/2} }{2}                              second approach 
//'
//' 
//' - xmin xmax - INPUT: pointers to arrays of double containing respectively the lower and upper
//'               boundaries of the different variables. For a given variable, only the values that are
//'               implied by the type of constraints, defined as in *type, are actually inspected. 
//' 
//' - f - function that calculates the objective function at a specified point x. Its
//'       specification is
//'       void (*f) (const size_t n, const double *x,void *fparams,double *fval)
//'       where
//'       n       - INPUT: the number of variables 
//'       x       - INPUT:the point at which the function is required 
//'       fparams - INPUT: pointer to a structure containing parameters required by the
//'                 function. If no external parameter are required it can be set to NULL. 


                    //'       fval    - OUTPUT: the value of the objective function at the current point x. 
//' 
//' - df - function that calculates the gradient of the objective function at a specified
//'        point x. Its specification is
//'        void (*df) (const size_t n, const double *x,void *fparams,double *grad)
//'        where
//'        n       - INPUT: the number of variables 
//'        x       - INPUT:the point at which the function is required 
//'        fparams - INPUT: pointer to a structure containing parameters required by the
//'                  function. If no external parameter are required it can be set to NULL. 
//'        grad    - OUTPUT: the values of the gradient of the objective function at the
//'                  current point x are stored in grad[0],...,grad[n-1]. 
//' 
//' - fdf - fdf calculates the value and the gradient of the objective function at a
//'         specified point x. Its specification is
//'         void (*fdf) (const size_t n, const double *x,void *fparams,double *fval,double *grad)
//'         where
//'         n       - INPUT: the number of variables 
//'         x       - INPUT:the point at which the function is required 
//'         fparams - INPUT: pointer to a structure containing parameters required by the
//'                   function. If no external parameter are required it can be set to NULL. 
//'         fval    - OUTPUT: the value of the objective function at the current point x. 
//'         grad    - OUTPUT: the values of the gradient of the objective function at the
//'                   current point x are stored in grad[0],...,grad[n-1]. 
//' 
//' - fparams - pointer to a structure containing parameters required by the function. If no
//'             external parameter are required it can be set to NULL. 
//' 
//' - oparams - structure of the type "multiminparams" containing the optimization
//'             parameters. The members are
//' 
//'             double step_size - size of the first trial step 
//'             double tol       - accuracy of the line minimization 
//'             unsigned maxiter - maximum number of iterations 
//'             double epsabs    - accuracy of the minimization 
//'             double maxsize   - final size of the simplex 
//'             unsigned method  - method to use. Possible values are:
//' 
//'              - Fletcher-Reeves conjugate gradient
//'              - Polak-Ribiere conjugate gradient
//'              - Vector Broyden-Fletcher-Goldfarb-Shanno method
//'              - Steepest descent algorithm
//'              - Nelder-Mead simplex
//'              - Vector Broyden-Fletcher-Goldfarb-Shanno ver. 2
//' 
//'             unsigned verbosity - if greater then 0 print info on intermediate steps 

void multimin(
               Rcpp::NumericVector data
              ,size_t n
	      ,std::vector<double> x
              ,double *fun
              ,Rcpp::IntegerVector type 
	      ,Rcpp::NumericVector xmin
	      ,Rcpp::NumericVector xmax
              ,void (*f)    (Rcpp::NumericVector, const size_t, std::vector<double>, void *, double *)
              ,void (* df)  (Rcpp::NumericVector, const size_t, std::vector<double>, void *, std::vector<double>)
              ,void (* fdf) (Rcpp::NumericVector, const size_t, std::vector<double>, void *, double *, std::vector<double>)
              ,void *fparams
              ,const struct multimin_params oparams
              ){


  size_t i=0;

  // choose algorithm for the optimization process
  struct multimin_algorithm multimin_alg  =
    choose_algorithm(oparams.method);
  
  gsl_vector *y  = gsl_vector_alloc (n);

  // selects the algoritm for the optimization procedure
  //const gsl_multimin_fdfminimizer_type *Tfdf = multimin_alg->Tfdf;
  //const gsl_multimin_fminimizer_type *Tf = multimin_alg->Tf;
  //const char *Tname = multimin_alg->Tname;
  
  
  /* --- OUPUT ---------------------------------- */
    Rprintf("MULTIMIN START\n");
    Rprintf("#    method                         %s\n", multimin_alg.Tname);

    if(oparams.method<4 || oparams.method==5){

      Rprintf("#    initial step size              %g\n", oparams.step_size);
      Rprintf("#    line minimization tolerance    %g\n", oparams.tol);
      Rprintf("#    maximum number of iterations   %u\n", oparams.maxiter);
      Rprintf("#    precision                      %g\n", oparams.epsabs);

    }
    else{

      Rprintf("#    maximum number of iterations   %u\n",oparams.maxiter);
      Rprintf("#    maximum simplex size           %g\n",oparams.maxsize);

    }
  /* -------------------------------------------- */

  // compute values of y for initial condition
  do_initial_data_transform( x
                            ,n
                            ,y
                            ,type
                            ,xmin
                            ,xmax
			    );


    {
      double res;
      Rprintf("Objective function initial value\n");
      f(data, n, x, fparams, &res);
      Rprintf("#    f=%.2e\n",res);
    }
  /* -------------------------------------------- */


  if(oparams.method<4 || oparams.method==5){/* methods with derivatives */

    unsigned iter=0;
    int status1,status2;
    
    // generic parameters struct
    struct g_params gparams;
    
    // generic function to be minimized
    gsl_multimin_function_fdf GdG;
    
    // solver to be used
    gsl_multimin_fdfminimizer *s = gsl_multimin_fdfminimizer_alloc (multimin_alg.Tfdf,n);

    /* set the parameters of the new function */
    gparams.data    = data;
    gparams.n       = n;
    gparams.type    = type;
    gparams.xmin    = xmin;
    gparams.xmax    = xmax;
    gparams.f       = f;
    gparams.df      = df;
    gparams.fdf     = fdf;
    gparams.fparams = fparams;
    
    /* set the function to solve */
    GdG.f=g;
    GdG.df=dg;
    GdG.fdf=gdg;
    GdG.n=n;
    GdG.params=(void *) &gparams;

  
    /* initialize minimizer */
    status1=gsl_multimin_fdfminimizer_set(s                        // algorithm used
                                          ,&GdG                    // function to be minimized
                                          ,y                       // initial point
                                          ,oparams.step_size       // step size
                                          ,oparams.tol             // tolerance of line minimization (error proxy)
                                          );  

    if(status1){
        Rcpp::stop("#ERROR: %s\n",gsl_strerror (status1));
    }

    /* +++++++++++++++++++++++++++++++++++++++++++++++ */
    if(oparams.verbosity>2){
       Rprintf("#    - start minimization \n");
    }
    /* +++++++++++++++++++++++++++++++++++++++++++++++ */
   
    do
      {

        if( ++iter > oparams.maxiter) break;
        
        status1 = gsl_multimin_fdfminimizer_iterate (s);

        /* +++++++++++++++++++++++++++++++++++++++++++++++ */
        if(oparams.verbosity>2){
          Rprintf("#     [%d]",iter);
          Rprintf(" g=%+12.6e  y=( ",s->f);
          for(i=0;i<n;i++)
            Rprintf("%+12.6e ",gsl_vector_get(s->x,i));
          Rprintf(") dg=( ");
          for(i=0;i<n;i++)
            Rprintf("%+12.6e  ",gsl_vector_get(s->gradient,i));
          Rprintf(") |dg|=%12.6e ",gsl_blas_dnrm2 (s->gradient));
          Rprintf("|dx|=%12.6e\n",gsl_blas_dnrm2 (s->dx));
        }
        /* +++++++++++++++++++++++++++++++++++++++++++++++ */


        if(status1 && status1 != GSL_ENOPROG){
          Rprintf("#WARNING: %s\n", gsl_strerror (status1));
          break;
        }

        /* status2 = gsl_multimin_test_gradient (s->gradient,oparams.epsabs);  */
        status2 = gsl_multimin_test_gradient (gsl_multimin_fdfminimizer_gradient(s),oparams.epsabs); 

        if(status1 == GSL_ENOPROG && status2==GSL_CONTINUE){
          Rprintf("#    status: %s\n",gsl_strerror (status1));
          break;
        }
        
     }
    while (status2 == GSL_CONTINUE);

    gsl_vector_memcpy (y,s->x);
    *fun=s->f;
    gsl_multimin_fdfminimizer_free (s);

    /* +++++++++++++++++++++++++++++++++++++++++++++++ */
    if(oparams.verbosity>2){
      Rprintf("#    - end minimization\n");
      Rprintf("#    iterations %u\n",iter-1);
    }
    /* +++++++++++++++++++++++++++++++++++++++++++++++ */

  }
  else{ /* methods without derivatives */

    unsigned iter=0;
    int status1,status2;
    double size;
    gsl_vector *ss = gsl_vector_alloc (n);
    struct g_params gparams;
    gsl_multimin_function G;
    gsl_multimin_fminimizer *s=gsl_multimin_fminimizer_alloc (multimin_alg.Tf,n);

    /* set the parameters of the new function */
    gparams.n       = n;
    gparams.type    = type;
    gparams.xmin    = xmin;
    gparams.xmax    = xmax;
    gparams.f       = f;
    gparams.fparams = fparams;

    /* set the function to solve */
    G.f=g;
    G.n=n;
    G.params=(void *) &gparams;

    /* Initial vertex size vector */
    gsl_vector_set_all (ss,oparams.step_size+oparams.maxsize);

    /* --- OUPUT ---------------------------------- */
    if(oparams.verbosity>0){
      size_t i;
      Rprintf("#    initial simplex sizes\n");
      Rprintf("#    ");
      for(i=0;i<n;i++)
        Rprintf(" %g", gsl_vector_get(ss,i));
      Rprintf("\n");
    }
    /* -------------------------------------------- */

    /* Initialize minimizer */ 
    status1=gsl_multimin_fminimizer_set(s,&G,y,ss);

    if(status1)
      {
        Rprintf("#ERROR: %s\n",gsl_strerror (status1));
        Rcpp::stop("Error on the minimization process.");
      }

    /* +++++++++++++++++++++++++++++++++++++++++++++++ */
    if(oparams.verbosity>2)
      Rprintf("#    - start minimization \n");
    /* +++++++++++++++++++++++++++++++++++++++++++++++ */

    do
      {

        if( ++iter > oparams.maxiter) break;

        status1 = gsl_multimin_fminimizer_iterate(s);

        size = gsl_multimin_fminimizer_size (s);

        /* +++++++++++++++++++++++++++++++++++++++++++++++ */
        if(oparams.verbosity>2){
          Rprintf("#    g=%g y=( ",s->fval);
          for(i=0;i<n;i++)
            Rprintf("%g ",gsl_vector_get(s->x,i));
          Rprintf(") ");
          Rprintf(" simplex size=%g ",size);
          Rprintf("\n");
        }
        /* +++++++++++++++++++++++++++++++++++++++++++++++ */


        if(status1 && status1 != GSL_ENOPROG){
          Rprintf("#WARNING: %s\n", gsl_strerror (status1));
          break;
        }

        status2=gsl_multimin_test_size (size,oparams.maxsize);

        if(status1 == GSL_ENOPROG && status2==GSL_CONTINUE){
          Rprintf("#    status: %s\n",gsl_strerror (status1));
          break;
        }

      }
    while (status2 == GSL_CONTINUE);

    gsl_vector_memcpy (y, s->x);
    *fun=s->fval;
    gsl_multimin_fminimizer_free (s);

    /* +++++++++++++++++++++++++++++++++++++++++++++++ */
    if(oparams.verbosity>2){
      Rprintf("#    - end minimization\n");
      Rprintf("#    iterations %u\n",iter-1);
    }
    /* +++++++++++++++++++++++++++++++++++++++++++++++ */

  }

  /* compute values of x */
  do_data_transformation(x, n, y, type, xmin, xmax);


  /* --- OUPUT ---------------------------------- */
  if(oparams.verbosity>0){
    for(i=0;i<n;i++)
      Rprintf("#    %e -> x[%zd]=%e\n",gsl_vector_get(y,i),i,x[i]);
    Rprintf("#--- MULTIMIN END --- \n");
  }
  /* -------------------------------------------- */


  gsl_vector_free (y);
  
}
