#include <Rcpp.h>
#include <assert.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector simpleSmoothing(NumericVector series, double alpha, double startValue) {
    if (alpha > 1.0 || alpha < 0.0)
        throw std::range_error("Inadmissible value");

    int length = series.length();
    double previous = startValue;

    NumericVector result = NumericVector::create(length);

    for (int i = 0; i < length; ++i) {
        result[i] = alpha * series[i] + (1.0 - alpha) * previous;
        previous = result[i];
    }

    return result;
}
