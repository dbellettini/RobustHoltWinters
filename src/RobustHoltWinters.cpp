/*
 *  R robust Holt-Winters extension
 *  Copyright (C) 2015 Davide Bellettini
 *
 *  based on R: A Computer Language for Statistical Data Analysis
 *  Copyright (C) 2003-7  The R Development Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.gnu.org/licenses/gpl-3.0.txt
 */

#include <Rcpp.h>
#include <assert.h>
#include <cmath>

// TODO: use these constant everywhere
// requires rstudio to update bindings
const double k = 2.0, ck = 2.52;

using namespace Rcpp;

double sgn(double x)
{
    if (x < 0) {
        return -1.0;
    }

    return 1.0;
}

/**
 * saturation function, returns x in the interval [-k, +k]
 * sign(x) * k otherwhise
 */
// [[Rcpp::export]]
double psi (double x)
{
    if (std::abs(x) > k) {
        return sgn(x) * k;
    }

    return x;
}

/*
 * res is the residual (e.g. x - xhat) thus the difference
 * between actual value and 1 step ahead prediction
 */
bool shouldSaturate(double residual, double sigma)
{
    return (residual / sigma) > k || (residual / sigma) < -k;
}

/*
 * rho is a bounded loss function that is used to replace
 * pow(x, 2) to get robust
 *
 * returns double [0, ck]
 */
double rho (double x)
{
    if (x > k || x < -k) {
        return ck;
    }

    return ck * (1 - pow(1 - pow(x / k, 2), 3));
}

/**
 * sigma estimates the scale of error and it's needed for the
 * saturation, and takes sigma up to time t-1
 */
double updatesigma (double delta, double rt, double sigma)
{
    return sqrt(
            delta * rho(rt / sigma) * pow(sigma, 2)
                + (1 - delta) * pow(sigma, 2)
    );
}

// [[Rcpp::export]]
List RobustHoltWintersCpp(
    NumericVector x,
    const double alpha,
    const double beta,
    const double gamma,
    int startTime,
    int frequency,
    double levelInitial,
    double trendInitial,
    NumericVector seasonInitial,
    double sigma,
    Function median
) {
    const int xLength = x.length();
    const int smoothedLength = xLength - startTime + 1;
    const int seasonLength = frequency;

    // sigma smoothing parameter (see the function updatesigma)
    double delta = 0.2;

    NumericVector level(smoothedLength + 1);
    NumericVector trend(smoothedLength + 1);
    NumericVector season(smoothedLength + seasonLength);

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
        const int previousIndex = currentIndex - 1;
        const int currentSeasonalIndex = currentIndex + seasonLength - 1;
        const double xAtT = x[t];

        /* forecast *for* period i */
        // xhat is the 1 step ahead prediction
        const double seasonalPrevious = season[currentSeasonalIndex - seasonLength];
        const double xhat = level[previousIndex] + trend[previousIndex] + seasonalPrevious;

        /* Sum of Squared Errors */
        double residual = xAtT - xhat;
        double xFiteredAtT = xAtT;

        if (shouldSaturate(residual, sigma)) {
            residual = sigma * psi(residual / sigma);
            xFiteredAtT = xhat + residual;
        }

        SSE += pow(residual, 2);

        /* estimate of level *in* period t */
        level[currentIndex] = alpha * (xFiteredAtT - seasonalPrevious)
            + (1 - alpha) * (level[previousIndex] + trend[previousIndex]);

        /* estimate of trend *in* period t */
        trend[currentIndex] = beta * (level[currentIndex] - level[previousIndex])
            + (1 - beta)  * trend[previousIndex];

        /* estimate of seasonal component *in* period t */
        season[currentSeasonalIndex] = gamma * (xFiteredAtT - level[currentIndex])
            + (1 - gamma) * seasonalPrevious;

        sigma = updatesigma(delta, residual, sigma);
    }

    List output;

    output["SSE"] = SSE;
    output["level"] = level;
    output["trend"] = trend;
    output["seasonal"] = season;

    return output;
}
