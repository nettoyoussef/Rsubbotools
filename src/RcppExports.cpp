// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppGSL.h>
#include <Rcpp.h>

using namespace Rcpp;

// varcovar
RcppGSL::Matrix varcovar(const Rcpp::NumericVector par, const size_t N, const size_t dim);
RcppExport SEXP _Rsubbotools_varcovar(SEXP parSEXP, SEXP NSEXP, SEXP dimSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector >::type par(parSEXP);
    Rcpp::traits::input_parameter< const size_t >::type N(NSEXP);
    Rcpp::traits::input_parameter< const size_t >::type dim(dimSEXP);
    rcpp_result_gen = Rcpp::wrap(varcovar(par, N, dim));
    return rcpp_result_gen;
END_RCPP
}
// sortRcpp
void sortRcpp(Rcpp::NumericVector x);
RcppExport SEXP _Rsubbotools_sortRcpp(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    sortRcpp(x);
    return R_NilValue;
END_RCPP
}
// subbofit
Rcpp::List subbofit(Rcpp::NumericVector data, int verb, int method, int interv_step, int output, Rcpp::Nullable<Rcpp::NumericVector> provided_m_, Rcpp::NumericVector par, Rcpp::NumericVector g_opt_par, Rcpp::NumericVector itv_opt_par);
RcppExport SEXP _Rsubbotools_subbofit(SEXP dataSEXP, SEXP verbSEXP, SEXP methodSEXP, SEXP interv_stepSEXP, SEXP outputSEXP, SEXP provided_m_SEXP, SEXP parSEXP, SEXP g_opt_parSEXP, SEXP itv_opt_parSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type verb(verbSEXP);
    Rcpp::traits::input_parameter< int >::type method(methodSEXP);
    Rcpp::traits::input_parameter< int >::type interv_step(interv_stepSEXP);
    Rcpp::traits::input_parameter< int >::type output(outputSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type provided_m_(provided_m_SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type par(parSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type g_opt_par(g_opt_parSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type itv_opt_par(itv_opt_parSEXP);
    rcpp_result_gen = Rcpp::wrap(subbofit(data, verb, method, interv_step, output, provided_m_, par, g_opt_par, itv_opt_par));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_Rsubbotools_varcovar", (DL_FUNC) &_Rsubbotools_varcovar, 3},
    {"_Rsubbotools_sortRcpp", (DL_FUNC) &_Rsubbotools_sortRcpp, 1},
    {"_Rsubbotools_subbofit", (DL_FUNC) &_Rsubbotools_subbofit, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_Rsubbotools(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
