#include <Rcpp.h>
#include <assert.h>
#include <cmath>

using namespace Rcpp;

// [[Rcpp::export]]
List SupervisedHoltWintersCpp(
    NumericVector x,
    IntegerVector outlier,
    const double alpha,
    const double beta,
    const double gamma,
    int startTime,
    int frequency,
    double levelInitial,
    double trendInitial,
    NumericVector seasonInitial
) {
    const int xLength = x.length();
    const int smoothedLength = xLength - startTime + 1;
    const int seasonLength = frequency;

    NumericVector level(smoothedLength + 1);
    NumericVector trend(smoothedLength + 1);
    NumericVector season(smoothedLength + seasonLength);
    NumericVector residuals(smoothedLength);

    /* copy start values to the beginning of the vectors */
    level[0] = levelInitial;
    trend[0] = trendInitial;

    for (int i = 0; i < seasonLength; ++i) {
        season[i] = seasonInitial[i];
    }

    double SSE = 0.0;

    for (int t = startTime - 1; t < xLength; t++) {
        /* indices for period i */
        const int currentIndex = t - startTime + 2;
        const int currentResidualIndex = currentIndex - 1;
        const int previousIndex = currentIndex - 1;
        const int currentSeasonalIndex = currentIndex + seasonLength - 1;
        const double xAtT = x[t];

        /* forecast *for* period i */
        // xhat is the 1 step ahead prediction
        const double seasonalPrevious =
            season[currentSeasonalIndex - seasonLength];
        const double xhat =
            level[previousIndex] + trend[previousIndex] + seasonalPrevious;

        /* Sum of Squared Errors */
        double residual = residuals[currentResidualIndex] = xAtT - xhat;
        double xFiteredAtT = xAtT;

        if (outlier[t]) {
            residual = 0;
            xFiteredAtT = xhat;
        }

        SSE += pow(residual, 2);

        /* estimate of level *in* period t */
        level[currentIndex] = alpha * (xFiteredAtT - seasonalPrevious)
            + (1 - alpha) * (level[previousIndex] + trend[previousIndex]);

        /* estimate of trend *in* period t */
        trend[currentIndex] =
            beta * (level[currentIndex] - level[previousIndex])
            + (1 - beta)  * trend[previousIndex];

        /* estimate of seasonal component *in* period t */
        season[currentSeasonalIndex] =
            gamma * (xFiteredAtT - level[currentIndex])
            + (1 - gamma) * seasonalPrevious;
    }

    List output;

    output["SSE"] = SSE;
    output["level"] = level;
    output["trend"] = trend;
    output["seasonal"] = season;

    return output;
}
