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

      int n = data.size();

      Rprintf("#>>> [%+.3e:%+.3e]", data[i], data[i+1]);
      for(int tmp=0; tmp < n;tmp++){
        Rprintf(" par[%d]= %e", tmp, par[i]);
      }
      Rprintf(" ll= %e\n", *fmin);
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

/* long options */
/* ------------ */

struct option gb_long_options[] = {
  {"version", no_argument,       NULL,  0 },
  {"help", no_argument,       NULL,  'h' },
  {0,         0,                 0,  0 }
};

int gb_option_index = 0;
