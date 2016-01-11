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
double psi (double x, double k)
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
bool shouldSaturate(double res, double sigma)
{
    return (res / sigma) > k || (res / sigma) < -k;
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
    int seasonalType,
    int frequency,
    double levelInitial,
    double trendInitial,
    NumericVector seasonInitial,
    double sigma,
    double k
) {
    int xl = x.length();
    int len = xl - startTime + 1;
    double delta = 0.2;

    double res = 0, xhat = 0, stmp = 0;
    int i, i0, s0;

    double SSE = 0.0;

    NumericVector level(len + 1);
    NumericVector trend(len + 1);
    NumericVector season(len + frequency);

    /* copy start values to the beginning of the vectors */
    level[0] = levelInitial;
    trend[0] = trendInitial;

    for (i = 0; i < frequency; ++i) {
        season[i] = seasonInitial[i];
    }

    for (i = startTime - 1; i < xl; i++) {
        /* indices for period i */
        i0 = i - startTime + 2;
        s0 = i0 + frequency - 1;

        /* forecast *for* period i */
        xhat = level[i0 - 1] + trend[i0 - 1];
        stmp = season[s0 - frequency];

        if (seasonalType == 1)
            xhat += stmp;
        else
            xhat *= stmp;

        /* Sum of Squared Errors */
        res   = x[i] - xhat;

        if (shouldSaturate(res, sigma)) {
            /* fprintf(stderr, "res: %f, sigma: %f, k: %f\n", res, sigma, k); */
            res = sigma * psi(res / sigma, k);
            x[i] = xhat + res;
        }

        SSE += res * res;

        /* estimate of level *in* period i */
        if (seasonalType == 1)
            level[i0] = alpha       * (x[i] - stmp)
                + (1 - alpha) * (level[i0 - 1] + trend[i0 - 1]);
        else
            level[i0] = alpha       * (x[i] / stmp)
                + (1 - alpha) * (level[i0 - 1] + trend[i0 - 1]);

        /* estimate of trend *in* period i */
        trend[i0] = beta        * (level[i0] - level[i0 - 1])
            + (1 - beta)  * trend[i0 - 1];

        /* estimate of seasonal component *in* period i */
        if (seasonalType == 1)
            season[s0] = gamma * (x[i] - level[i0])
                + (1 - gamma) * stmp;
        else
            season[s0] = gamma * (x[i] / level[i0])
                + (1 - gamma) * stmp;

        sigma = updatesigma(delta, res, sigma);
    }

    List output;

    output["SSE"] = SSE;
    output["level"] = level;
    output["trend"] = trend;
    output["seasonal"] = season;

    return output;
}
