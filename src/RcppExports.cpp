// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// RobustHoltWintersCpp
List RobustHoltWintersCpp(NumericVector x, NumericVector filtered, const double alpha, const double beta, const double gamma, int startTime, int seasonalType, int frequency, bool doTrend, bool doSeasonal, double levelInitial, double trendInitial, NumericVector seasonInitial);
RcppExport SEXP RobustHoltWinters_RobustHoltWintersCpp(SEXP xSEXP, SEXP filteredSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP gammaSEXP, SEXP startTimeSEXP, SEXP seasonalTypeSEXP, SEXP frequencySEXP, SEXP doTrendSEXP, SEXP doSeasonalSEXP, SEXP levelInitialSEXP, SEXP trendInitialSEXP, SEXP seasonInitialSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type filtered(filteredSEXP );
        Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP );
        Rcpp::traits::input_parameter< const double >::type beta(betaSEXP );
        Rcpp::traits::input_parameter< const double >::type gamma(gammaSEXP );
        Rcpp::traits::input_parameter< int >::type startTime(startTimeSEXP );
        Rcpp::traits::input_parameter< int >::type seasonalType(seasonalTypeSEXP );
        Rcpp::traits::input_parameter< int >::type frequency(frequencySEXP );
        Rcpp::traits::input_parameter< bool >::type doTrend(doTrendSEXP );
        Rcpp::traits::input_parameter< bool >::type doSeasonal(doSeasonalSEXP );
        Rcpp::traits::input_parameter< double >::type levelInitial(levelInitialSEXP );
        Rcpp::traits::input_parameter< double >::type trendInitial(trendInitialSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type seasonInitial(seasonInitialSEXP );
        List __result = RobustHoltWintersCpp(x, filtered, alpha, beta, gamma, startTime, seasonalType, frequency, doTrend, doSeasonal, levelInitial, trendInitial, seasonInitial);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
