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
  Code explanation

  Below there is an overview of the GSL Multimin Suite.
  This suite as a front end for the user has several shortcomings:

  - it requires different minimization procedures for methods with
    derivatives and methods without derivatives, e.g., changing function calls
    and names.

  - it is unable to restrict the parameters range, which
    requires that variables be transformed previously to being minimized.
    This is cumbersome, as it forces the user to change the objective function
    for each change in the range of acceptable parameters.

  Both these restrictions make it very redudant to change across methods
  with derivatives from methods without or restricting the range of parameters
  when one is using the same objective function.

  The multmin mask, which occurs on this routine, overcomes this problem.
  In my opinion, this should be the default in the gsl library, and I am not
  sure why it was not made so.

  Basically, what is done is:
  - parameters transformation - this transformation suitable limits the values
    that the parameters can have before being passed to the minimizer.
  - masks for the objective functions f, df, fdf (see below), so that all
    functions have the same arguments and same style.
  - a single routine that decides by checking the optimization method which
    initialization routine should be used.

  The multimin call

  void multimin(
                Rcpp::NumericVector data
                ,size_t n
                ,Rcpp::NumericVector x
                ,double *fun
                ,Rcpp::IntegerVector type
                ,Rcpp::NumericVector xmin
                ,Rcpp::NumericVector xmax
                ,void (*f)    (Rcpp::NumericVector, const size_t, Rcpp::NumericVector, void *, double *)
                ,void (* df)  (Rcpp::NumericVector, const size_t, Rcpp::NumericVector, void *, Rcpp::NumericVector)
                ,void (* fdf) (Rcpp::NumericVector, const size_t, Rcpp::NumericVector, void *, double *, Rcpp::NumericVector)
                ,void *fparams
                ,const struct multimin_params oparams
                ,int verb
                ){


  The mask requires that the objective functions created by the user have
  the following format:

    void (*f)    (Rcpp::NumericVector data, const size_t n, Rcpp::NumericVector x, void *params, double *f)
    void (* df)  (Rcpp::NumericVector data, const size_t n, Rcpp::NumericVector x, void *params, Rcpp::NumericVector df)
    void (* fdf) (Rcpp::NumericVector data, const size_t n, Rcpp::NumericVector x, void *params, double *f, Rcpp::NumericVector df)

    - data   - the sample of points for the function to be minimized over
    - n      -  size of the sample
    - x      - vector with current value of the parameters
    - params - other params that are not hardcoded on the objective function, but are not optimized
    - *f     - the objective function value, which is updated by reference
    - *df    - the gradient of the objective function, which is updated by reference

  Other parameters for the multimin function:
    - Rcpp::IntegerVector type - vector containing the ranges
    - Rcpp::NumericVector xmin - vector containing the minimum value for the parameters
    - Rcpp::NumericVector xmax - vector containing the maximum value for the parameters
    - void *fparams            - other parameters passed for the function
    - const struct multimin_params oparams - multimin optimization parameters
    - int verb                 - verbosity level

  The function works by creating masks of the basic functions f, df, fdf which
  suits the format expected by GSL, respectively g, dg, gdg. Those functions
  receive the original functions as a argument on the struct gparams.

  This struct is basically composed of:

    gparams.data    = data;
    gparams.n       = n;
    gparams.type    = type;
    gparams.xmin    = xmin;
    gparams.xmax    = xmax;
    gparams.f       = f;
    gparams.df      = df;
    gparams.fdf     = fdf;
    gparams.fparams = fparams;

   This mask works by first applying the reverse transformation to the data,
   given by the function do_initial_data_transform




  GSL Multimin minimization suite - an overview

  Note: Methods that minimize functions using derivatives use functions with the
  fdf prefix, while methods without derivatives use only an f, e.g.
  fminimizer vs fdfminimizer.
  Traditionally, the GSL suite requires the following code to work:

   - gsl_multimin_fdfminimizer *pointer / gsl_multimin_fminimizer *pointer
      the minimizer struct, allocated as a pointer;

   - const gsl_multimin_fdfminimizer_type *T / const gsl_multimin_fminimizer_type *T
      the minimizer type allocated as a pointer. It can be one of the following:

       0: [conjugate_fr    ] Fletcher-Reeves conjugate gradient
       1: [conjugate_pr    ] Polak-Ribiere conjugate gradient
       2: [vector_bfgs     ] Vector Broyden-Fletcher-Goldfarb-Shanno method
       3: [steepest_descent] Steepest descent algorithm
       4: [nmsimplex       ] Nelder-Mead simplex
       5: [vector_bfgs2    ] Vector Broyden-Fletcher-Goldfarb-Shanno method ver. 2
       6: [nmsimplex2      ] Simplex algorithm of Nelder and Mead ver. 2
       7: [nmsimplex2rand  ] Simplex algorithm of Nelder and Mead: random initialization

       For the codes, just add the value between brackets to the prefix gsl_multimin_fminimizer_*.
       Options 0-3 and 5 require derivatives, while 4, 6 and 7 don't.
       Methods with derivatives use the gsl_multimin_fdfminimizer_type, while
       methods without use the gsl_multimin_fminimizer_type.

       This vector needs to be initiatilized by the following functions
       - gsl_multimin_fdfminimizer_alloc (T, n) - where T is the above type,
         already allocated with a suitable algorithm, and n is the size of the
         parameters' vector x.

       Example:
       const gsl_multimin_fdfminimizer_type *T;
       T = gsl_multimin_fminimizer_vector_bfgs2;
       gsl_multimin_fdfminimizer *s;
       s = gsl_multimin_fdfminimizer_alloc (T, n);

   - gsl_multimin_function_fdf my_func - the minimizer function struct which must
     respect the following struct:

      my_func.f     = double (* f) (const gsl_vector * x, void * params);
      my_func.df    = void (* df) (const gsl_vector * x, void * params, gsl_vector * g);
      my_func.fdf   = void (* fdf) (const gsl_vector * x, void * params, double * f, gsl_vector * g);
      my_func.n     = size_t n;
      my_func.params= void * params;

      - f      - the function that the minimizer aims to optimize, must return a value or
                 GSL_NAN.
      - df     - this function must store the gradient for argument x on vector g
                 or return an error if this is not possible.
      - fdf    - this function updates by reference the values of f and g above.
                 It is usually the preferred method since it seems to be faster.
      - n      - the number of parameters to be optimized, i.e., the size of vector x.
      - params - other parameters that are not optimized but are passed to the
                 function.

   - gsl_multimin_function my_func - for methods without derivatives the
     structure is:

      my_func.f     = double (* f) (const gsl_vector * x, void * params);
      my_func.n     = size_t n;
      my_func.params= void * params;

   - gsl_vector *x - the pointer which holds the initial guess, and that will be updated
     by the algorithm iteration.

       This vector needs to be initiatilized by the following functions
       - gsl_vector_alloc (n) - allocates in the memory a vector of size n
       - gsl_vector_set(x, pos, val) - fills the vector in the position pos with the value val
       Example:
       gsl_vector *x;
       x = gsl_vector_alloc (2);
       gsl_vector_set (x, 0, 5.0);
       gsl_vector_set (x, 1, 7.0);

   - gsl_multimin_fdfminimizer_set ( gsl_multimin_fdfminimizer *s
                                    ,gsl_multimin_function_fdf *fdf
                                    ,const gsl_vector *x
                                    ,double step_size
                                    ,double tol
                                    )
     gsl_multimin_fminimizer_set( gsl_multimin_fminimizer * s
                                 ,gsl_multimin_function * f
                                 ,const gsl_vector * x
                                 ,const gsl_vector * step_size
                                 )

     These functions initialize the minimizer s with the function f/fdf, the initial
     guess x, the initial step size step_size and, when available, the tolerance tol.

   - gsl_multimin_fminimizer_iterate / gsl_multimin_fdfminimizer_iterate
     performs one iteration over the minimizer s, after it was initialized above.

     The struct s has the following elements:

     For the fdfminimizer
     s->type     - minimizer type;
     s->x        - vector with positions of the current guess;
     s->dx       - vector with last step increments over current guess
     s->f        - objective function value
     s->gradient - vector with gradient values
     s->state    - error codes
     s->fdf      - function used?

     For the fminimizer
     s->type     - minimizer type;
     s->x        - vector with positions of the current guess;
     s->f        - used in this context?
     s->fval     - objective function value
     s->size     - minimizer specific characteristic size
     s->state    - error codes



     ### A full example for a method without derivative, from the docs: ###
     int main (void){
       size_t iter = 0;
       int status;

       const gsl_multimin_fdfminimizer_type *T;
       gsl_multimin_fdfminimizer *s;

       // Position of the minimum (1,2), scale factors
          10,20, height 30. //
       double par[5] = { 1.0, 2.0, 10.0, 20.0, 30.0 };

       gsl_vector *x;
       gsl_multimin_function_fdf my_func;

       my_func.n = 2;
       my_func.f = my_f;
       my_func.df = my_df;
       my_func.fdf = my_fdf;
       my_func.params = par;

       // Starting point, x = (5,7) //
       x = gsl_vector_alloc (2);
       gsl_vector_set (x, 0, 5.0);
       gsl_vector_set (x, 1, 7.0);

       T = gsl_multimin_fdfminimizer_conjugate_fr;
       s = gsl_multimin_fdfminimizer_alloc (T, 2);

       gsl_multimin_fdfminimizer_set (s, &my_func, x, 0.01, 1e-4);

       do
         {
           iter++;
           status = gsl_multimin_fdfminimizer_iterate (s);

           if (status)
             break;

           status = gsl_multimin_test_gradient (s->gradient, 1e-3);

           if (status == GSL_SUCCESS)
             printf ("Minimum found at:\n");

           printf ("%5d %.5f %.5f %10.5f\n", iter,
                   gsl_vector_get (s->x, 0),
                   gsl_vector_get (s->x, 1),
                   s->f);

         }
       while (status == GSL_CONTINUE && iter < 100);

       gsl_multimin_fdfminimizer_free (s);
       gsl_vector_free (x);

       return 0;
     }

     ### A full example for a method with derivative, from the docs: ###

     int main(void){
       double par[5] = {1.0, 2.0, 10.0, 20.0, 30.0};

       const gsl_multimin_fminimizer_type *T =
         gsl_multimin_fminimizer_nmsimplex2;
       gsl_multimin_fminimizer *s = NULL;
       gsl_vector *ss, *x;
       gsl_multimin_function minex_func;

       size_t iter = 0;
       int status;
       double size;

       // Starting point //
       x = gsl_vector_alloc (2);
       gsl_vector_set (x, 0, 5.0);
       gsl_vector_set (x, 1, 7.0);

       // Set initial step sizes to 1 //
       ss = gsl_vector_alloc (2);
       gsl_vector_set_all (ss, 1.0);

       // Initialize method and iterate //
       minex_func.n = 2;
       minex_func.f = my_f;
       minex_func.params = par;

       s = gsl_multimin_fminimizer_alloc (T, 2);
       gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

       do
         {
           iter++;
           status = gsl_multimin_fminimizer_iterate(s);

           if (status)
             break;

           size = gsl_multimin_fminimizer_size (s);
           status = gsl_multimin_test_size (size, 1e-2);

           if (status == GSL_SUCCESS)
             {
               printf ("converged to minimum at\n");
             }

           printf ("%5d %10.3e %10.3e f() = %7.3f size = %.3f\n",
                   iter,
                   gsl_vector_get (s->x, 0),
                   gsl_vector_get (s->x, 1),
                   s->fval, size);
         }
       while (status == GSL_CONTINUE && iter < 100);

       gsl_vector_free(x);
       gsl_vector_free(ss);
       gsl_multimin_fminimizer_free (s);

       return status;
     }
*/

#include "multimin.h"

// compute values of y for initial condition
RcppGSL::vector<double> do_initial_data_transform(
                               Rcpp::NumericVector x
                               ,const size_t n
                               ,Rcpp::IntegerVector type
                               ,Rcpp::NumericVector xmin
                               ,Rcpp::NumericVector xmax
                               ,int verb
                               ){


  RcppGSL::vector<double> y(n);
  size_t i;
  double dtmp1;

  if(verb > 1){
    Rprintf("INITIAL DATA TRANSFORMATION\n");
    Rprintf("#    - variables initial value and boundaries\n");
  }


  for(i=0;i<n;i++){
    if(type[i] == NA_INTEGER){
      y[i] = x[i];
        if(verb > 1){
          Rprintf("#    x[%d]=%.2f (-inf,+inf)   trans 0 -> %.2f\n",(int) i,x[i],gsl_vector_get(y,i));
        }
    }
    else
      switch(type[i]){
      case 0:/* (-inf,+inf) */
        y[i] = x[i];
        if(verb > 1){
          Rprintf("#    x[%d]=%.2f (-inf,+inf)   trans 0 -> %.2f\n"
                  ,(int) i,x[i],gsl_vector_get(y,i));
        }
        break;
      case 1:/* [a,+inf) */
        y[i] = sqrt( x[i]-xmin[i] );
        if(verb > 1){
          Rprintf("#    x[%d]=%.2f [%.3g,+inf)   trans 1 -> %.2f\n"
                  ,(int) i,x[i],xmin[i],gsl_vector_get(y,i));
        }
        break;
      case 2:/* (-inf,a] */
        y[i] = sqrt( xmax[i]-x[i] );
        if(verb > 1){
          Rprintf("#    x[%d]=%.2f (-inf,%.3f]        trans 2 -> %.2f\n"
                  ,(int) i,x[i],xmax[i],gsl_vector_get(y,i));
        }
        break;
      case 3:/* [a,b] */
        dtmp1 = (xmax[i]>xmin[i]?
                 (2.*x[i]-xmax[i]-xmin[i])/(xmax[i]-xmin[i]) : 0);
        /*       dtmp1 = (2.*x[i]-xmax[i]-xmin[i])/(xmax[i]-xmin[i]); */
        y[i] = asin( dtmp1 );
        if(verb > 1){
          Rprintf("#    x[%d]=%.2f [%.3f,%.3f] trans 3 -> %.2f\n"
                  ,(int) i,x[i],xmin[i],xmax[i],gsl_vector_get(y,i));
        }
        break;
      case 4:/* (a,+inf) */
        y[i] = log( x[i]-xmin[i] );
        if(verb > 1){
          Rprintf("#    x[%d]=%.2f (%.3f,+inf)  trans 4 -> %.2f\n"
                  ,(int) i,x[i],xmin[i],gsl_vector_get(y,i));
        }
        break;
      case 5:/* (-inf,a) */
        y[i] = log( xmax[i]-x[i] );
        if(verb > 1){
          Rprintf("#    x[%d]=%.2f (-inf,%.3f)  trans 5 -> %.2f\n"
                  ,(int) i,x[i],xmax[i],gsl_vector_get(y,i));
        }
        break;
      case 6:/* (a,b) */
        dtmp1 = (2.*x[i]-xmax[i]-xmin[i])/(xmax[i]-xmin[i]);
        y[i] = gsl_atanh ( dtmp1 );
        if(verb > 1){
          Rprintf("#    x[%d]=%.2f (%.3f,%.3f) trans 6 -> %.2f\n"
                  ,(int) i,x[i],xmin[i],xmax[i],gsl_vector_get(y,i));
        }
        break;
      case 7:/* (a,b) second approach */
        dtmp1 = (2.*x[i]-xmax[i]-xmin[i])/(xmax[i]-xmin[i]);
        y[i] =  dtmp1/sqrt(1-dtmp1*dtmp1);
        if(verb > 1){
          Rprintf("#    x[%d]=%.2f (%.3f,%.3f) trans 7 -> %.2f\n"
                  ,(int) i,x[i],xmin[i],xmax[i],gsl_vector_get(y,i));
        }
        break;
      case 8:/* (a,+inf) second approach */
        dtmp1 = x[i]-xmin[i];
        y[i] =  dtmp1-1./(4.*dtmp1);
        if(verb > 1){
          Rprintf("#    x[%d]=%.2f (%.3f,+inf)  trans 8 -> %.2f\n"
                  ,(int) i,x[i],xmin[i],gsl_vector_get(y,i));
        }
        break;
      case 9:/* (-inf,a) second approach */
        dtmp1 = xmax[i]-x[i];
        y[i] =  1./(4.*dtmp1)-dtmp1;
        if(verb > 1){
          Rprintf("#    x[%d]=%.2f (-inf,%.3f)  trans 9 -> %.2f\n"
                  ,(int) i,x[i],xmax[i],gsl_vector_get(y,i));
        }
        break;
      }
  }

  if(verb > 1){
    Rprintf("END OF DATA TRANSFORMATION\n");
  }

  return y;
}


void do_data_transformation(
                            Rcpp::NumericVector x
                            ,const size_t n
                            ,Rcpp::NumericVector y
                            ,Rcpp::IntegerVector type
                            ,Rcpp::NumericVector xmin
                            ,Rcpp::NumericVector xmax
                            ){


  size_t i;
  double dtmp1;

  for(i=0;i<n;i++){
    if(type[i] == NA_INTEGER){
      x[i]= y[i];
    }
    else
      switch(type[i]){
      case 0:/* (-inf,+inf) */
        x[i]= y[i];
        break;
      case 1:/* [a,+inf) */
        x[i]= xmin[i]+y[i]*y[i];
        break;
      case 2:/* (-inf,a] */
        x[i]= xmax[i]-y[i]*y[i];
        break;
      case 3:/* [a,b] */
        dtmp1 = sin( y[i] );
        x[i]= .5*(xmin[i]*(1-dtmp1) +xmax[i]*(1+dtmp1));
        break;
      case 4:/* (a,+inf) */
        dtmp1 = exp( y[i] );
        x[i]= xmin[i]+dtmp1;
        break;
      case 5:/* (-inf,a) */
        dtmp1 = -exp( y[i] );
        x[i]= xmax[i]+dtmp1;
        break;
      case 6:/* (a,b) */
        dtmp1 = tanh( y[i] );
        x[i]= .5*(xmin[i]*(1-dtmp1) +xmax[i]*(1+dtmp1));
        break;
      case 7:/* (a,b) second approach */
        dtmp1 = y[i]/sqrt(1.+y[i]*y[i]);
        x[i]= .5*(xmin[i]*(1-dtmp1) +xmax[i]*(1+dtmp1));
        break;
      case 8:/* (a,+inf) second approach */
        dtmp1 = sqrt(1.+y[i]*y[i]);
        x[i]= xmin[i] + .5*(y[i]+dtmp1);
        break;
      case 9:/* (-inf,a) second approach */
        dtmp1 = sqrt(1.+y[i]*y[i]);
        x[i]= xmax[i] + .5*(y[i]-dtmp1);
        break;
      }
  }
}


void do_data_transformation_df(
                               Rcpp::NumericVector dx
                               ,const size_t n
                               ,Rcpp::NumericVector y
                               ,Rcpp::IntegerVector type
                               ,Rcpp::NumericVector xmin
                               ,Rcpp::NumericVector xmax
                               ){

  //Rcpp::NumericVector x
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
        dx[i]= 2.*y[i];
        break;
      case 2:/* (-inf,a] */
        dx[i]= -2.*y[i];
        break;
      case 3:/* [a,b] */
        dx[i]= .5*(xmax[i]-xmin[i])*cos(y[i]);
        break;
      case 4:/* (a,+inf) */
        dx[i]= exp( y[i] );
        break;
      case 5:/* (-inf,a) */
        dx[i]= -exp( y[i] );
        break;
      case 6:/* (a,b) */
        dtmp1= cosh( y[i] );
        dx[i]= .5*(xmax[i]-xmin[i])/(dtmp1*dtmp1);
        break;
      case 7:/* (a,b) second approach */
        dtmp1= (1.+y[i]*y[i])*sqrt(1.+y[i]*y[i]) ;
        dx[i]= .5*(xmax[i]-xmin[i])/dtmp1;
        break;
      case 8:/* (a,+inf) second approach */
        dtmp1= y[i];
        dtmp2= sqrt(1.+dtmp1*dtmp1);
        dx[i]= .5*(dtmp1+dtmp2)/dtmp2;
        break;
      case 9:/* (-inf,a) second approach */
        dtmp1 = y[i];
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

  // the function is forced to return a value
  double res=GSL_NAN;
  struct g_params *p= (struct g_params *) gparams;

  /* dereference useful stuff */
  const size_t n           = p->n;
  Rcpp::IntegerVector type = p->type;
  Rcpp::NumericVector xmin = p->xmin;
  Rcpp::NumericVector xmax = p->xmax;

  //double *x = (double *) multimin_alloc(sizeof(double)*n);
  Rcpp::NumericVector x(n);
  //std::vector<double> x(n);

  // converts y from a gls_vector to a Rcpp vector
  Rcpp::NumericVector z = Rcpp::NumericVector(y->data, y->data + y->size);
  //Rcpp::Rcout << "Print z:" << std::endl;
  //Rcpp::Rcout << z << std::endl;

  // variable transformation

  // type is used to select the boundary conditions for
  // the values of the different variables/coefficents.
  // the process uses a data transformation to convert a
  // constrained problem in a unconstrained one.

  // Notice: type is also a vector, so for each
  // variable a different transformation can be selected

  // data transf. for x
  do_data_transformation(x, n, z, type, xmin, xmax);

  //Rcpp::Rcout << "Print x:" << std::endl;
  //Rcpp::Rcout << x << std::endl;

  //updates the objective function
  // res - variable with this value, passed by reference
  p->f(p->data, n, x, p->fparams,&res);

  //Rcpp::Rcout << "Print res:" << std::endl;
  //Rcpp::Rcout << res << std::endl;

  //free (x);
  return(res);

}



// this function must allocate in a pointer the
// value of the gradient of f evaluated at point
// y_1, ..., y_n
static void dg(const gsl_vector *y, void *gparams, gsl_vector *dg){

  struct g_params *p= (struct g_params *) gparams;

  // converts y to a Rcpp vector
  Rcpp::NumericVector z = Rcpp::NumericVector(y->data, y->data + y->size);

  /* dereference useful stuff */
  const size_t n = p->n;
  Rcpp::IntegerVector type = p->type;
  Rcpp::NumericVector xmin = p->xmin;
  Rcpp::NumericVector xmax = p->xmax;

  // x is the guess after the variable transformation
  Rcpp::NumericVector x(n);
  //std::vector<double> x(n);
  //double *x = (double *) multimin_alloc(sizeof(double)*n);

  // dx is the gradient after the variable transformation
  Rcpp::NumericVector dx(n);
  //std::vector<double> dx(n);
  //double *dx = (double *) multimin_alloc(sizeof(double)*n);

  Rcpp::NumericVector df(n);
  //std::vector<double> df(n);
  //double *df = (double *) multimin_alloc(sizeof(double)*n);




  // variable transformation

  // type is used to select the boundary conditions for
  // the values of the different variables/coefficents.
  // the process uses a data transformation to convert a
  // constrained problem in a unconstrained one.

  // Notice: type is also a vector, so for each
  // variable a different transformation can be selected

  // data transf. for x
  do_data_transformation(x, n, z, type, xmin, xmax);

  // data transf. for dx
  do_data_transformation_df(dx, n, z, type, xmin, xmax);

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

  // converts y to a Rcpp vector
  Rcpp::NumericVector z = Rcpp::NumericVector(y->data, y->data + y->size);

  /* dereference useful stuff */
  const size_t n = p->n;
  Rcpp::IntegerVector type = p->type;
  Rcpp::NumericVector xmin = p->xmin;
  Rcpp::NumericVector xmax = p->xmax;


  // x is the guess after the variable transformation
  Rcpp::NumericVector x(n);
  //std::vector<double> x(n);
  //double *x = (double *) multimin_alloc(sizeof(double)*n);

  // dx is the gradient after the variable transformation
  Rcpp::NumericVector dx(n);
  //std::vector<double> dx(n);
  //double *dx = (double *) multimin_alloc(sizeof(double)*n);

  Rcpp::NumericVector df(n);
  //std::vector<double> df(n);
  //double *df = (double *) multimin_alloc(sizeof(double)*n);

  // variable transformation

  // type is used to select the boundary conditions for
  // the values of the different variables/coefficents.
  // the process uses a data transformation to convert a
  // constrained problem in a unconstrained one.

  // Notice: type is also a vector, so for each
  // variable a different transformation can be selected

  // data transf. for x
  do_data_transformation(x, n, z, type, xmin, xmax);

  // data transf. for dx
  do_data_transformation_df(dx, n, z, type, xmin, xmax);

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
struct multimin_algorithm choose_algorithm(unsigned int method
                                           ,int verb
                                           ){


  if(verb > 1){
    Rprintf("Choosing algorithm:\n");
  }

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

  if(verb > 1){
    Rprintf("Algorithm chosen: %s\n", multimin_alg.Tname);
  }

  return multimin_alg;

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



//' multimin
//' Interface for the GSL Optmization routine
//'
//'
//'
//' @param size_t n - dimension of the problem, number of independent parameters of
//' the objective function
//' @param double *x - pointer used as input for the array of n values x[0],...
//' x[n-1] containing the initial estimate of the minimum point. The pointer is
//' updated with the final estimation of the minimum position.
//' @param double *fun -
//' @param unsigned *type - a pointer to an array of integer
//' type[1],...,type[n-1] describing the boundary conditions for the different
//' variables. The problem is solved as an unconstrained one on a suitably
//' transformed variable y. Possible values range from integers 0-9. See details below.
//' @param double *xmin, *xmax - pointers to arrays of double containing
//' respectively the lower and upper boundaries of the different n parameters.
//' For a given parameter, only the values that are implied by the type of
//' constraints, defined as in *type, are actually inspected.
//' @param void (*f) - objective function. See details below.
//' @param void (* df) - gradient of the objective function, for the methods that
//' use it. See details below.
//' @param void (* fdf) - combination of the objective function and the gradient.
//' Some methods may be faster to be calculated when intermediate steps are calculated
//' for the function and gradient at the same function. GSL checks this, and use
//' the faster version.
//' @param void *fparams - special parameters that are not optimized and are
//' requested by the objective function. If not necessary, use NULL.
//' @param const struct multimin_params oparams - struct containing
//'
//' --------------------------------------------------------------------
//' type
//'
//' Possible values are:
//'   Interval:                                       Transformation:
//'   0 unconstrained                                 x=y
//'   1 semi-closed right half line [ xmin,+infty )   x=xmin+y^2
//'   2 semi-closed left  half line ( -infty,xmax ]   x=xmax-y^2
//'   3 closed interval              [ xmin,xmax ]    x=SS+SD*sin(y)
//'   4 open right half line        ( xmin,+infty )   x=xmin+exp(y)
//'   5 open left  half line        ( -infty,xmax )   x=xmax-exp(y)
//'   6 open interval                ( xmin,xmax )    x=SS+SD*tanh(y)
//'
//' where SS=.5(xmin+xmax) SD=.5(xmax-xmin)
//'
//' There are also other UNSUPPORTED transformations used in various test
//'   7 open interval                ( xmin,xmax )    x=SS+SD*(1+y/sqrt(1+y^2))
//'   8 open right half line        ( xmin,+infty )   x=xmin+.5*(y+sqrt(1+y^2))
//'   9 open left  half line        ( -infty,xmax )   x=xmax+.5*(y-sqrt(1+y^2))
//' --------------------------------------------------------------------
//' f
//'
//' f calculates the objective function at a specified point x. Its
//' specification is
//'
//' void (*f) (const size_t n, const double *x,void *fparams,double *fval)
//'
//'       n
//'       INPUT: the number of variables
//'
//'       x
//'       INPUT:the point at which the function is required
//'
//'       fparams
//'       pointer to a structure containing parameters required by the
//'       function. If no external parameter are required it can be set to
//'       NULL.
//'
//'       fval
//'       OUTPUT: the value of the objective function at the current point
//'       x.
//'
//' --------------------------------------------------------------------
//' df
//'
//' df calculates the gradient of the objective function at a specified
//' point x. Its specification is
//'
//' void (*df) (const size_t n, const double *x,void *fparams,double *grad)
//'
//'       n
//'       INPUT: the number of variables
//'
//'       x
//'       INPUT:the point at which the function is required
//'
//'       fparams
//'       pointer to a structure containing parameters required by the
//'       function. If no external parameter are required it can be set to
//'       NULL.
//'
//'       grad
//'       OUTPUT: the values of the gradient of the objective function at
//'       the current point x are stored in grad[0],...,grad[n-1].
//'
//' --------------------------------------------------------------------
//' fdf
//'
//' fdf calculates the value and the gradient of the objective function at
//' a specified point x. Its specification is
//'
//' void (*fdf) (const size_t n, const double *x,void *fparams,double *fval,double *grad)
//'
//'       n
//'       INPUT: the number of variables
//'
//'       x
//'       INPUT:the point at which the function is required
//'
//'       fparams
//'       pointer to a structure containing parameters required by the
//'       function. If no external parameter are required it can be set to
//'       NULL.
//'
//'       fval
//'       OUTPUT: the value of the objective function at the current point
//'       x.
//'
//'       grad
//'       OUTPUT: the values of the gradient of the objective function at
//'       the current point x are stored in grad[0],...,grad[n-1].
//'
//' --------------------------------------------------------------------
//' fparams
//'
//' pointer to a structure containing parameters required by the
//' function. If no external parameter are required it can be set to NULL.
//'
//' --------------------------------------------------------------------
//'
//' oparams
//'
//' structure of the type "multimin_params" containing the optimization
//' parameters. The members are
//'
//'       double step_size
//'           size of the first trial step
//'
//'       double tol
//'           accuracy of the line minimization
//'
//'       unsigned maxiter
//'           maximum number of iterations
//'
//'       double epsabs
//'           accuracy of the minimization
//'
//'       double maxsize;
//'           final size of the simplex
//'
//'       unsigned method
//'           method to use. Possible values are:
//'
//'           0: Fletcher-Reeves conjugate gradient
//'           1: Polak-Ribiere conjugate gradient
//'           2: Vector Broyden-Fletcher-Goldfarb-Shanno method
//'           3: Steepest descent algorithm
//'           4: Nelder-Mead simplex
//'           5: Vector Broyden-Fletcher-Goldfarb-Shanno method ver. 2
//'           6: Simplex algorithm of Nelder and Mead ver. 2
//'           7: Simplex algorithm of Nelder and Mead: random initialization
//'
//'       unsigned verbosity
//'           if greater then 0 print info on intermediate steps

void multimin(
              Rcpp::NumericVector data
              ,size_t n
              ,Rcpp::NumericVector x
              ,double *fun
              ,Rcpp::IntegerVector type
              ,Rcpp::NumericVector xmin
              ,Rcpp::NumericVector xmax
              ,void (*f)    (Rcpp::NumericVector, const size_t, Rcpp::NumericVector, void *, double *)
              ,void (* df)  (Rcpp::NumericVector, const size_t, Rcpp::NumericVector, void *, Rcpp::NumericVector)
              ,void (* fdf) (Rcpp::NumericVector, const size_t, Rcpp::NumericVector, void *, double *, Rcpp::NumericVector)
              ,void *fparams
              ,const struct multimin_params oparams
              ,int verb
              ){


  size_t i=0;

  // choose algorithm for the optimization process
  struct multimin_algorithm multimin_alg  =
    choose_algorithm(oparams.method, verb);

  RcppGSL::vector<double> y(n);


  /* --- OUPUT ---------------------------------- */
  if(verb > 1){

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

  }
  /* -------------------------------------------- */
  // objective function initial value
  {
    double res;

    if(verb > 1){
      Rprintf("# Objective function initial value\n");
    }
    f(data, n, x, fparams, &res);
    if(verb > 1){
      Rprintf("# f = %.2e\n", res);
    }
  }

  // make transformation on the x ranges compute values of y for initial condition
  y  = do_initial_data_transform(x, n, type, xmin, xmax, verb);

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
    GdG.f      = g;
    GdG.df     = dg;
    GdG.fdf    = gdg;
    GdG.n      = n;
    GdG.params = (void *) &gparams;

    /* initialize minimizer */
    status1=gsl_multimin_fdfminimizer_set(s                        // minimizer
                                          ,&GdG                    // function to be minimized
                                          ,y                       // initial point
                                          ,oparams.step_size       // step size
                                          ,oparams.tol             // tolerance of line minimization (error proxy)
                                          );

    if(status1){
      Rcpp::stop("#ERROR: %s\n",gsl_strerror (status1));
    }

    /* +++++++++++++++++++++++++++++++++++++++++++++++ */
    if(verb > 1){
      Rprintf("#    - start minimization \n");
    }
    //}
    /* +++++++++++++++++++++++++++++++++++++++++++++++ */

    do
      {

        if( ++iter > oparams.maxiter) break;

        status1 = gsl_multimin_fdfminimizer_iterate (s);

        /* +++++++++++++++++++++++++++++++++++++++++++++++ */
        if(verb > 1){

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
    if(verb > 1){
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
    gsl_vector_set_all (ss, oparams.step_size+oparams.maxsize);

    /* --- OUPUT ---------------------------------- */
    if(verb > 1){
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
    if(verb > 1){
      Rprintf("#    - start minimization \n");
    }
    /* +++++++++++++++++++++++++++++++++++++++++++++++ */

    do
      {

        if( ++iter > oparams.maxiter) break;

        status1 = gsl_multimin_fminimizer_iterate(s);
        size = gsl_multimin_fminimizer_size (s);

        /* +++++++++++++++++++++++++++++++++++++++++++++++ */
        if(verb > 1){
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
    if(verb > 1){
      Rprintf("#    - end minimization\n");
      Rprintf("#    iterations %u\n",iter-1);
    }
    /* +++++++++++++++++++++++++++++++++++++++++++++++ */

    }

  /* compute values of x */
  Rcpp::NumericVector z = Rcpp::NumericVector(y->data, y->data + y->size);
  do_data_transformation(x, n, z, type, xmin, xmax);


  /* --- OUPUT ---------------------------------- */
  if(verb > 1){
    for(i=0;i<n;i++)
      Rprintf("#    %e -> x[%zd]=%e\n",gsl_vector_get(y,i),i,x[i]);
    Rprintf("#--- MULTIMIN END --- \n");
  }
  /* -------------------------------------------- */

  // gsl_vector_free (y);

}
