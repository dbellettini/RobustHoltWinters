#include <Rcpp.h>
#include <assert.h>
#include <iostream>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector simpleSmoothing(NumericVector series, const double alpha, const double startValue) {
    if (alpha > 1.0 || alpha < 0.0)
        throw std::range_error("Inadmissible value");

    const int length = series.length();
    NumericVector result(length);

    double previous = startValue;

    for (int i = 0; i < length; ++i) {
        result[i] = alpha * series[i] + (1.0 - alpha) * previous;
        previous = result[i];
    }

    return result;
}
