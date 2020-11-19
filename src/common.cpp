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


