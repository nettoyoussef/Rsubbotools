/*insert GNU extensions*/
#define _GNU_SOURCE
/*in particular, use of NAN extension*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <assert.h>
#include <getopt.h>
//#include <config.h> // to check

// R messages
#include<Rcpp.h>

/* used by <errno.h> */
extern int errno;

/* GSL ----------- */
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>

/* --------------- */

// #include "config.h"

void *my_alloc(size_t);

void *my_realloc(void *,size_t);

int sort_by_value(const void *, const void *);

int load(double **,unsigned *,FILE *);

void moment(double *, double *, double *, double *);

void subbotools_header(char const *,FILE *);


/* C++ ----------- */
#include<RcppGSL.h>
// [[Rcpp::depends(RcppGSL)]]

#include<Rcpp.h>
/*- -------------- */

// conversion functions
double * RcppNum_to_double(Rcpp::NumericVector x);

// utility functions

// sort a vector
void sortRcpp(Rcpp::NumericVector x);

// copy a vector by value
void Rcppdeepcopy(Rcpp::NumericVector x_orig, Rcpp::NumericVector x_dest);

// functions for optimizing the Subbotin Family fit
#include "multimin.h"

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
                          ,void (*f)    (Rcpp::NumericVector, const size_t, Rcpp::NumericVector, void *, double *)
                          ,void (* df)  (Rcpp::NumericVector, const size_t, Rcpp::NumericVector, void *, Rcpp::NumericVector)
                          ,void (* fdf) (Rcpp::NumericVector, const size_t, Rcpp::NumericVector, void *, double *, Rcpp::NumericVector)
                          ,int verb
                          );

Rcpp::List global_optim(
                        Rcpp::NumericVector data
                        ,double fmin
                        ,struct multimin_params global_oparams
                        ,Rcpp::NumericVector par
                        ,unsigned n_param
                        ,unsigned m_position
                        ,void (*f)    (Rcpp::NumericVector, const size_t, Rcpp::NumericVector, void *, double *)
                        ,void (* df)  (Rcpp::NumericVector, const size_t, Rcpp::NumericVector, void *, Rcpp::NumericVector)
                        ,void (* fdf) (Rcpp::NumericVector, const size_t, Rcpp::NumericVector, void *, double *, Rcpp::NumericVector)
                        ,Rcpp::Nullable<Rcpp::NumericVector> provided_m_
                        ,int verb
                        );


/* functions for the information matrix of asymmetric subbotin  */
double B0(double);
double B1(double);
double B2(double);
double dB0dx(double);
double dB0dx2(double);

/* long options */
extern struct option gb_long_options[];
extern int gb_option_index;
