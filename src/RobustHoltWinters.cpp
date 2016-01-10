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
#include <iostream>
#include <cmath>

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
double psi (double y, double k)
{
    if (y > k) {
        return k;
    }

    if (y < -k) {
        return -k;
    }

    return y;
}

bool shouldSaturate(double res, double sigma, double k)
{
    return (res / sigma) > k || (res / sigma) < -k;
}

double rho (double x)
{
    const double k = 2, ck = 2.52;

    if (x > k || x < -k) {
        return ck;
    }

    return ck * (1 - pow(1 - pow(x / k, 2), 3));
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
    bool doTrend,
    bool doSeasonal,
    double levelInitial,
    double trendInitial,
    NumericVector seasonInitial,
    double sigma,
    double k
) {
    int xl = x.length();
    int len = xl - startTime + 1;

    double res = 0, xhat = 0, stmp = 0;
    int i, i0, s0;

    double SSE = 0.0;

    NumericVector level(len + 1);
    NumericVector trend(len + 1);
    NumericVector season(len + frequency);

    /* copy start values to the beginning of the vectors */
    level[0] = levelInitial;
    if (doTrend) {
        trend[0] = trendInitial;
    }

    if (doSeasonal) {
        for (i = 0; i < frequency; ++i) {
            season[i] = seasonInitial[i];
        }
    }

    for (i = startTime - 1; i < xl; i++) {
        /* indices for period i */
        i0 = i - startTime + 2;
        s0 = i0 + frequency - 1;

        /* forecast *for* period i */
        xhat = level[i0 - 1] + (doTrend ? trend[i0 - 1] : 0);
        stmp = doSeasonal ? season[s0 - frequency] : (seasonalType != 1);

        if (seasonalType == 1)
            xhat += stmp;
        else
            xhat *= stmp;

        /* Sum of Squared Errors */
        res   = x[i] - xhat;

        if (shouldSaturate(res, sigma, k)) {
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
        if (doTrend)
            trend[i0] = beta        * (level[i0] - level[i0 - 1])
                + (1 - beta)  * trend[i0 - 1];

        /* estimate of seasonal component *in* period i */
        if (doSeasonal) {
            if (seasonalType == 1)
                season[s0] = gamma * (x[i] - level[i0])
                    + (1 - gamma) * stmp;
            else
                season[s0] = gamma * (x[i] / level[i0])
                    + (1 - gamma) * stmp;
        }
    }

    List output;

    output["SSE"] = SSE;
    output["level"] = level;

    if (doTrend) {
        output["trend"] = trend;
    }

    if (doSeasonal) {
        output["seasonal"] = season;
    }

    return output;
}
