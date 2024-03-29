// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppGSL.h>
#include <Rcpp.h>

using namespace Rcpp;

// rgamma
Rcpp::NumericVector rgamma(unsigned n, const double a, const double b);
RcppExport SEXP _Rsubbotools_rgamma(SEXP nSEXP, SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned >::type n(nSEXP);
    Rcpp::traits::input_parameter< const double >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(rgamma(n, a, b));
    return rcpp_result_gen;
END_RCPP
}
// rlaplace
Rcpp::NumericVector rlaplace(unsigned n, const double m, const double a);
RcppExport SEXP _Rsubbotools_rlaplace(SEXP nSEXP, SEXP mSEXP, SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned >::type n(nSEXP);
    Rcpp::traits::input_parameter< const double >::type m(mSEXP);
    Rcpp::traits::input_parameter< const double >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(rlaplace(n, m, a));
    return rcpp_result_gen;
END_RCPP
}
// ralaplace
Rcpp::NumericVector ralaplace(unsigned n, const double m, const double al, const double ar);
RcppExport SEXP _Rsubbotools_ralaplace(SEXP nSEXP, SEXP mSEXP, SEXP alSEXP, SEXP arSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned >::type n(nSEXP);
    Rcpp::traits::input_parameter< const double >::type m(mSEXP);
    Rcpp::traits::input_parameter< const double >::type al(alSEXP);
    Rcpp::traits::input_parameter< const double >::type ar(arSEXP);
    rcpp_result_gen = Rcpp::wrap(ralaplace(n, m, al, ar));
    return rcpp_result_gen;
END_RCPP
}
// rpower
Rcpp::NumericVector rpower(unsigned n, const double a, const double b);
RcppExport SEXP _Rsubbotools_rpower(SEXP nSEXP, SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned >::type n(nSEXP);
    Rcpp::traits::input_parameter< const double >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(rpower(n, a, b));
    return rcpp_result_gen;
END_RCPP
}
// rsubbo
Rcpp::NumericVector rsubbo(unsigned n, double m, double a, double b);
RcppExport SEXP _Rsubbotools_rsubbo(SEXP nSEXP, SEXP mSEXP, SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(rsubbo(n, m, a, b));
    return rcpp_result_gen;
END_RCPP
}
// rasubbo
Rcpp::NumericVector rasubbo(unsigned n, double m, double bl, double br, double al, double ar);
RcppExport SEXP _Rsubbotools_rasubbo(SEXP nSEXP, SEXP mSEXP, SEXP blSEXP, SEXP brSEXP, SEXP alSEXP, SEXP arSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type bl(blSEXP);
    Rcpp::traits::input_parameter< double >::type br(brSEXP);
    Rcpp::traits::input_parameter< double >::type al(alSEXP);
    Rcpp::traits::input_parameter< double >::type ar(arSEXP);
    rcpp_result_gen = Rcpp::wrap(rasubbo(n, m, bl, br, al, ar));
    return rcpp_result_gen;
END_RCPP
}
// alaplafit
Rcpp::List alaplafit(Rcpp::NumericVector data, int verb, unsigned interv_step, Rcpp::Nullable<Rcpp::NumericVector> provided_m_);
RcppExport SEXP _Rsubbotools_alaplafit(SEXP dataSEXP, SEXP verbSEXP, SEXP interv_stepSEXP, SEXP provided_m_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type verb(verbSEXP);
    Rcpp::traits::input_parameter< unsigned >::type interv_step(interv_stepSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type provided_m_(provided_m_SEXP);
    rcpp_result_gen = Rcpp::wrap(alaplafit(data, verb, interv_step, provided_m_));
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
// laplafit
Rcpp::List laplafit(Rcpp::NumericVector data, int verb, unsigned interv_step, Rcpp::Nullable<Rcpp::NumericVector> provided_m_);
RcppExport SEXP _Rsubbotools_laplafit(SEXP dataSEXP, SEXP verbSEXP, SEXP interv_stepSEXP, SEXP provided_m_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type verb(verbSEXP);
    Rcpp::traits::input_parameter< unsigned >::type interv_step(interv_stepSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type provided_m_(provided_m_SEXP);
    rcpp_result_gen = Rcpp::wrap(laplafit(data, verb, interv_step, provided_m_));
    return rcpp_result_gen;
END_RCPP
}
// sepfit
Rcpp::List sepfit(Rcpp::NumericVector data, int verb, Rcpp::NumericVector par, Rcpp::NumericVector g_opt_par);
RcppExport SEXP _Rsubbotools_sepfit(SEXP dataSEXP, SEXP verbSEXP, SEXP parSEXP, SEXP g_opt_parSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type verb(verbSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type par(parSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type g_opt_par(g_opt_parSEXP);
    rcpp_result_gen = Rcpp::wrap(sepfit(data, verb, par, g_opt_par));
    return rcpp_result_gen;
END_RCPP
}
// subboafish
Rcpp::List subboafish(unsigned size, double bl, double br, double m, double al, double ar, unsigned O_munknown);
RcppExport SEXP _Rsubbotools_subboafish(SEXP sizeSEXP, SEXP blSEXP, SEXP brSEXP, SEXP mSEXP, SEXP alSEXP, SEXP arSEXP, SEXP O_munknownSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< double >::type bl(blSEXP);
    Rcpp::traits::input_parameter< double >::type br(brSEXP);
    Rcpp::traits::input_parameter< double >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type al(alSEXP);
    Rcpp::traits::input_parameter< double >::type ar(arSEXP);
    Rcpp::traits::input_parameter< unsigned >::type O_munknown(O_munknownSEXP);
    rcpp_result_gen = Rcpp::wrap(subboafish(size, bl, br, m, al, ar, O_munknown));
    return rcpp_result_gen;
END_RCPP
}
// subboafit
Rcpp::List subboafit(Rcpp::NumericVector data, int verb, int method, int interv_step, Rcpp::Nullable<Rcpp::NumericVector> provided_m_, Rcpp::NumericVector par, Rcpp::NumericVector g_opt_par, Rcpp::NumericVector itv_opt_par);
RcppExport SEXP _Rsubbotools_subboafit(SEXP dataSEXP, SEXP verbSEXP, SEXP methodSEXP, SEXP interv_stepSEXP, SEXP provided_m_SEXP, SEXP parSEXP, SEXP g_opt_parSEXP, SEXP itv_opt_parSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type verb(verbSEXP);
    Rcpp::traits::input_parameter< int >::type method(methodSEXP);
    Rcpp::traits::input_parameter< int >::type interv_step(interv_stepSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type provided_m_(provided_m_SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type par(parSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type g_opt_par(g_opt_parSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type itv_opt_par(itv_opt_parSEXP);
    rcpp_result_gen = Rcpp::wrap(subboafit(data, verb, method, interv_step, provided_m_, par, g_opt_par, itv_opt_par));
    return rcpp_result_gen;
END_RCPP
}
// subbofish
Rcpp::List subbofish(unsigned size, double b, double m, double a, unsigned O_munknown);
RcppExport SEXP _Rsubbotools_subbofish(SEXP sizeSEXP, SEXP bSEXP, SEXP mSEXP, SEXP aSEXP, SEXP O_munknownSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< unsigned >::type O_munknown(O_munknownSEXP);
    rcpp_result_gen = Rcpp::wrap(subbofish(size, b, m, a, O_munknown));
    return rcpp_result_gen;
END_RCPP
}
// mm
double mm(const double std_over_aad, int verb);
RcppExport SEXP _Rsubbotools_mm(SEXP std_over_aadSEXP, SEXP verbSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type std_over_aad(std_over_aadSEXP);
    Rcpp::traits::input_parameter< int >::type verb(verbSEXP);
    rcpp_result_gen = Rcpp::wrap(mm(std_over_aad, verb));
    return rcpp_result_gen;
END_RCPP
}
// optim_method_moments
Rcpp::List optim_method_moments(Rcpp::NumericVector data, double fmin, Rcpp::Nullable<Rcpp::NumericVector> provided_m_, int verb);
RcppExport SEXP _Rsubbotools_optim_method_moments(SEXP dataSEXP, SEXP fminSEXP, SEXP provided_m_SEXP, SEXP verbSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type data(dataSEXP);
    Rcpp::traits::input_parameter< double >::type fmin(fminSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type provided_m_(provided_m_SEXP);
    Rcpp::traits::input_parameter< int >::type verb(verbSEXP);
    rcpp_result_gen = Rcpp::wrap(optim_method_moments(data, fmin, provided_m_, verb));
    return rcpp_result_gen;
END_RCPP
}
// subbofit
Rcpp::List subbofit(Rcpp::NumericVector data, int verb, int method, int interv_step, Rcpp::Nullable<Rcpp::NumericVector> provided_m_, Rcpp::NumericVector par, Rcpp::NumericVector g_opt_par, Rcpp::NumericVector itv_opt_par);
RcppExport SEXP _Rsubbotools_subbofit(SEXP dataSEXP, SEXP verbSEXP, SEXP methodSEXP, SEXP interv_stepSEXP, SEXP provided_m_SEXP, SEXP parSEXP, SEXP g_opt_parSEXP, SEXP itv_opt_parSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type verb(verbSEXP);
    Rcpp::traits::input_parameter< int >::type method(methodSEXP);
    Rcpp::traits::input_parameter< int >::type interv_step(interv_stepSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type provided_m_(provided_m_SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type par(parSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type g_opt_par(g_opt_parSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type itv_opt_par(itv_opt_parSEXP);
    rcpp_result_gen = Rcpp::wrap(subbofit(data, verb, method, interv_step, provided_m_, par, g_opt_par, itv_opt_par));
    return rcpp_result_gen;
END_RCPP
}
// subbolafit
Rcpp::List subbolafit(Rcpp::NumericVector data, int verb, int method, int interv_step, int output, Rcpp::Nullable<Rcpp::NumericVector> provided_m_, Rcpp::NumericVector par, Rcpp::NumericVector g_opt_par, Rcpp::NumericVector itv_opt_par);
RcppExport SEXP _Rsubbotools_subbolafit(SEXP dataSEXP, SEXP verbSEXP, SEXP methodSEXP, SEXP interv_stepSEXP, SEXP outputSEXP, SEXP provided_m_SEXP, SEXP parSEXP, SEXP g_opt_parSEXP, SEXP itv_opt_parSEXP) {
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
    rcpp_result_gen = Rcpp::wrap(subbolafit(data, verb, method, interv_step, output, provided_m_, par, g_opt_par, itv_opt_par));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_Rsubbotools_rgamma", (DL_FUNC) &_Rsubbotools_rgamma, 3},
    {"_Rsubbotools_rlaplace", (DL_FUNC) &_Rsubbotools_rlaplace, 3},
    {"_Rsubbotools_ralaplace", (DL_FUNC) &_Rsubbotools_ralaplace, 4},
    {"_Rsubbotools_rpower", (DL_FUNC) &_Rsubbotools_rpower, 3},
    {"_Rsubbotools_rsubbo", (DL_FUNC) &_Rsubbotools_rsubbo, 4},
    {"_Rsubbotools_rasubbo", (DL_FUNC) &_Rsubbotools_rasubbo, 6},
    {"_Rsubbotools_alaplafit", (DL_FUNC) &_Rsubbotools_alaplafit, 4},
    {"_Rsubbotools_sortRcpp", (DL_FUNC) &_Rsubbotools_sortRcpp, 1},
    {"_Rsubbotools_laplafit", (DL_FUNC) &_Rsubbotools_laplafit, 4},
    {"_Rsubbotools_sepfit", (DL_FUNC) &_Rsubbotools_sepfit, 4},
    {"_Rsubbotools_subboafish", (DL_FUNC) &_Rsubbotools_subboafish, 7},
    {"_Rsubbotools_subboafit", (DL_FUNC) &_Rsubbotools_subboafit, 8},
    {"_Rsubbotools_subbofish", (DL_FUNC) &_Rsubbotools_subbofish, 5},
    {"_Rsubbotools_mm", (DL_FUNC) &_Rsubbotools_mm, 2},
    {"_Rsubbotools_optim_method_moments", (DL_FUNC) &_Rsubbotools_optim_method_moments, 4},
    {"_Rsubbotools_subbofit", (DL_FUNC) &_Rsubbotools_subbofit, 8},
    {"_Rsubbotools_subbolafit", (DL_FUNC) &_Rsubbotools_subbolafit, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_Rsubbotools(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
