#  File src/library/stats/R/HoltWinters.R
#  Part of the R package, http://www.R-project.org
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

# Originally contributed by David Meyer

RobustHoltWinters <-
function (x, seasonal = c("additive", "multiplicative"))
{
    # starting values for optim
    optim.start = c(alpha = 0.3, beta = 0.1, gamma = 0.1)
    optim.control = list()
    start.periods = 2
    k = 2

    x <- as.ts(x)
    seasonal <- match.arg(seasonal)
    f <- frequency(x)

    ## initialization
    ## seasonal Holt-Winters
    start.time <- f + 1
    wind       <- start.periods * f

    ## decompose series
    st <- decompose(ts(x[1L:wind], start = start(x), frequency = f),
                    seasonal)

    ## level & intercept
    dat <- na.omit(st$trend)
    m   <- lm(dat ~ seq_along(dat))

    l.start <- as.vector(coef(m)[1L])
    b.start <- as.vector(coef(m)[2L])
    s.start <- st$figure

    ## Call to filtering loop
    lenx <- as.integer(length(x))
    if (is.na(lenx)) stop("invalid length(x)")

    len <- lenx - start.time + 1
    hw <- function (alpha, beta, gamma) {
        return (RobustHoltWintersCpp(
            x,
            as.double(alpha),
            as.double(beta),
            as.double(gamma),
            as.integer(start.time),
            as.integer(f),
            l.start,
            b.start,
            s.start,
            mad(x[1:start.time - 1]),
            k
        ))
    }

    ## --> optimize alpha, beta, and gamma
    error <- function (p) hw(p[1L], p[2L], p[3L])$SSE
    sol   <- optim(optim.start, error, method = "L-BFGS-B",
                   lower = c(0, 0, 0), upper = c(1, 1, 1),
                   control = optim.control)
    if(sol$convergence || any(sol$par < 0 | sol$par > 1)) {
        if (sol$convergence > 50) {
            warning(gettextf("optimization difficulties: %s",
                             sol$message), domain = NA)
        } else stop("optimization failure")
    }
    alpha <- sol$par[1L]
    beta  <- sol$par[2L]
    gamma <- sol$par[3L]

    ## get (final) results
    final.fit <- hw(alpha, beta, gamma)

    ## return fitted values and estimated coefficients along with parameters used
    fitted <- ts(cbind(xhat   = final.fit$level[-len-1],
                       level  = final.fit$level[-len-1],
                       trend  = if (!is.logical(beta) || beta)
                           final.fit$trend[-len-1],
                       season = if (!is.logical(gamma) || gamma)
                           final.fit$seasonal[1L:len]),
                 start = start(lag(x, k = 1 - start.time)),
                 frequency  = frequency(x)
                 )
    if (!is.logical(beta) || beta) fitted[,1] <- fitted[,1] + fitted[,"trend"]
    if (!is.logical(gamma) || gamma)
      fitted[,1] <- if (seasonal == "multiplicative")
        fitted[,1] * fitted[,"season"]
      else
        fitted[,1] + fitted[,"season"]
    structure(list(fitted    = fitted,
                   x         = x,
                   alpha     = alpha,
                   beta      = beta,
                   gamma     = gamma,
                   coefficients = c(a = final.fit$level[len + 1],
                                    b = if (!is.logical(beta) || beta) final.fit$trend[len + 1],
                                    s = if (!is.logical(gamma) || gamma) final.fit$seasonal[len + 1L:f]),
                   seasonal  = seasonal,
                   SSE       = final.fit$SSE,
                   call      = match.call()
                   ),
              class = "HoltWinters"
              )
}
