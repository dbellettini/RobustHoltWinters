// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// psi
double psi(double x);
RcppExport SEXP RobustHoltWinters_psi(SEXP xSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< double >::type x(xSEXP );
        double __result = psi(x);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// RobustHoltWintersCpp
List RobustHoltWintersCpp(NumericVector x, const double alpha, const double beta, const double gamma, int startTime, int frequency, double levelInitial, double trendInitial, NumericVector seasonInitial, double sigma, Function st);
RcppExport SEXP RobustHoltWinters_RobustHoltWintersCpp(SEXP xSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP gammaSEXP, SEXP startTimeSEXP, SEXP frequencySEXP, SEXP levelInitialSEXP, SEXP trendInitialSEXP, SEXP seasonInitialSEXP, SEXP sigmaSEXP, SEXP stSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP );
        Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP );
        Rcpp::traits::input_parameter< const double >::type beta(betaSEXP );
        Rcpp::traits::input_parameter< const double >::type gamma(gammaSEXP );
        Rcpp::traits::input_parameter< int >::type startTime(startTimeSEXP );
        Rcpp::traits::input_parameter< int >::type frequency(frequencySEXP );
        Rcpp::traits::input_parameter< double >::type levelInitial(levelInitialSEXP );
        Rcpp::traits::input_parameter< double >::type trendInitial(trendInitialSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type seasonInitial(seasonInitialSEXP );
        Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP );
        Rcpp::traits::input_parameter< Function >::type st(stSEXP );
        List __result = RobustHoltWintersCpp(x, alpha, beta, gamma, startTime, frequency, levelInitial, trendInitial, seasonInitial, sigma, st);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// SupervisedHoltWintersCpp
List SupervisedHoltWintersCpp(NumericVector x, IntegerVector outlier, const double alpha, const double beta, const double gamma, int startTime, int frequency, double levelInitial, double trendInitial, NumericVector seasonInitial);
RcppExport SEXP RobustHoltWinters_SupervisedHoltWintersCpp(SEXP xSEXP, SEXP outlierSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP gammaSEXP, SEXP startTimeSEXP, SEXP frequencySEXP, SEXP levelInitialSEXP, SEXP trendInitialSEXP, SEXP seasonInitialSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type outlier(outlierSEXP );
        Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP );
        Rcpp::traits::input_parameter< const double >::type beta(betaSEXP );
        Rcpp::traits::input_parameter< const double >::type gamma(gammaSEXP );
        Rcpp::traits::input_parameter< int >::type startTime(startTimeSEXP );
        Rcpp::traits::input_parameter< int >::type frequency(frequencySEXP );
        Rcpp::traits::input_parameter< double >::type levelInitial(levelInitialSEXP );
        Rcpp::traits::input_parameter< double >::type trendInitial(trendInitialSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type seasonInitial(seasonInitialSEXP );
        List __result = SupervisedHoltWintersCpp(x, outlier, alpha, beta, gamma, startTime, frequency, levelInitial, trendInitial, seasonInitial);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
