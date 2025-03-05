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
#include<R.h>
# include <Rmath.h>     // pgamma()
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
// multimin already includes structs, which contains data structures for
// all routines
//#include "structs.h"

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


// optimization functions for quantiles
double newton_c(newton_args x);
double steffensen_c(newton_args x);

/* functions for the information matrix of asymmetric subbotin  */
double B0(double);
double B1(double);
double B2(double);
double dB0dx(double);
double dB0dx2(double);

/* long options */
extern struct option gb_long_options[];
extern int gb_option_index;

// functions for the Laplace fit
double median(Rcpp::NumericVector data, size_t size);
double calculate_index(size_t size);

// functions for the quantile and RNG
double inc_lower_gamma(double b, double p);
double inv_inc_lower_gamma(double b, double p);
double inc_upper_gamma(double b, double p);
double inv_inc_upper_gamma(double b, double p);


// Templates

// returns sign of variable
// source
// https://stackoverflow.com/a/4609795/7233796
// templates must be on header files
// see
// https://stackoverflow.com/q/5612791/7233796
template <typename T> int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}


template <typename T>
    std::vector<std::size_t> rank_vector(const std::vector<T>& vector){

    // source:  https://stackoverflow.com/a/30685882/7233796

    // vector that is sorted, with same size as the original vector
    std::vector<std::size_t> indices(vector.size());

    // increases the vector by +1 starting on zero
    std::iota(indices.begin(), indices.end(), 0u);

    // sort with custom function
    // this function returns the original position of the entries
    // from the sorted vector
    std::sort(indices.begin(), indices.end(), [&](int lhs, int rhs) {
      // the sort criteria is the left element of the original vector
      // being smaller than the right one
      return vector[lhs] < vector[rhs];
    });

    // create another empty vector
    std::vector<std::size_t> res(vector.size());

    for (std::size_t i = 0; i < indices.size(); ++i) {
      res[indices[i]] = i;
    }

    //for (std::size_t i = 0; i < indices.size(); ++i) {
    //  Rprintf("indices[%d] =%d, res[%d] = %d\n", i, indices[i], i, res[i]);
    //}

    return res;
    //return indices;
  }
