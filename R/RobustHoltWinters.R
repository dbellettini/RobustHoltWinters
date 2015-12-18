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

RobustHoltWinters <-
function (x,

          # smoothing parameters
          alpha    = NULL, # level
          beta     = NULL, # trend
          gamma    = NULL, # seasonal component
          seasonal = c("additive", "multiplicative"),
          start.periods = 2,

          # starting values
          l.start  = NULL, # level
          b.start  = NULL, # trend
          s.start  = NULL, # seasonal components vector of length `period'

          # starting values for optim
          optim.start = c(alpha = 0.3, beta = 0.1, gamma = 0.1),
          optim.control = list()
          )
{
    prefilter <- identity

    x <- as.ts(x)
    prefiltered <- prefilter(x)
    seasonal <- match.arg(seasonal)
    f <- frequency(x)

    if(!is.null(alpha) && (alpha == 0))
        stop ("cannot fit models without level ('alpha' must not be 0 or FALSE).")
    if(!all(is.null(c(alpha, beta, gamma))) &&
        any(c(alpha, beta, gamma) < 0 || c(alpha, beta, gamma) > 1))
        stop ("'alpha', 'beta' and 'gamma' must be within the unit interval.")
    if((is.null(gamma) || gamma > 0)) {
        if (seasonal == "multiplicative" && any(x == 0))
            stop ("data must be non-zero for multiplicative Holt-Winters.")
        if (start.periods < 2)
            stop ("need at least 2 periods to compute seasonal start values.")
    }

    ## initialization
    if(!is.null(gamma) && is.logical(gamma) && !gamma) {
        ## non-seasonal Holt-Winters
        expsmooth <- !is.null(beta) && is.logical(beta) && !beta
        if(is.null(l.start))
            l.start <- if(expsmooth) x[1L] else x[2L]
        if(is.null(b.start))
            if(is.null(beta) || !is.logical(beta) || beta)
                b.start <- x[2L] - x[1L]
        start.time <- 3 - expsmooth
        s.start    <- 0
    } else {
        ## seasonal Holt-Winters
        start.time <- f + 1
        wind       <- start.periods * f

        ## decompose series
        st <- decompose(ts(x[1L:wind], start = start(x), frequency = f),
                        seasonal)

        ## level & intercept
        dat <- na.omit(st$trend)
        m   <- lm(dat ~ seq_along(dat))

        if (is.null(l.start)) l.start <- as.vector(coef(m)[1L])
        if (is.null(b.start)) b.start <- as.vector(coef(m)[2L])
        if (is.null(s.start)) s.start <- st$figure
    }

    ## Call to filtering loop
    len <- length(x) - start.time + 1

    hw <- function (alpha, beta, gamma) {
        return (RobustHoltWintersCpp(
            x,
            prefiltered,
            as.double(alpha),
            as.double(beta),
            as.double(gamma),
            as.integer(start.time),
            as.integer(! + (seasonal == "multiplicative")),
            as.integer(f),
            !is.logical(beta) || beta,
            !is.logical(gamma) || gamma,
            l.start,
            b.start,
            s.start
        ))
    }

    ## if alpha and/or beta and/or gamma are omitted, use optim to find the
    ## values minimizing the squared prediction error
    if (is.null(gamma)) {
        ## optimize gamma
        if (is.null(alpha)) {
            ## optimize alpha
            if (is.null(beta)) {
                ## optimize beta
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
            } else {
                ## !optimize beta
                ## --> optimize alpha and gamma
                error <- function (p) hw(p[1L], beta, p[2L])$SSE
                sol   <- optim(c(optim.start["alpha"], optim.start["gamma"]),
                               error, method = "L-BFGS-B",
                               lower = c(0, 0), upper = c(1, 1),
                               control = optim.control)
                if(sol$convergence || any(sol$par < 0 | sol$par > 1)) {
                    if (sol$convergence > 50) {
                        warning(gettextf("optimization difficulties: %s",
                                         sol$message), domain = NA)
                    } else stop("optimization failure")
                }
                alpha <- sol$par[1L]
                gamma <- sol$par[2L]
            }
        } else {
            ## !optimize alpha
            if (is.null(beta)) {
                ## optimize beta
                ## --> optimize beta and gamma
                error <- function (p) hw(alpha, p[1L], p[2L])$SSE
                sol   <- optim(c(optim.start["beta"], optim.start["gamma"]),
                               error, method = "L-BFGS-B",
                               lower = c(0, 0), upper = c(1, 1),
                               control = optim.control)
                if(sol$convergence || any(sol$par < 0 | sol$par > 1)) {
                    if (sol$convergence > 50) {
                        warning(gettextf("optimization difficulties: %s",
                                         sol$message), domain = NA)
                    } else stop("optimization failure")
                }
                beta  <- sol$par[1L]
                gamma <- sol$par[2L]
            } else {
                ## !optimize beta
                ## --> optimize gamma
                error <- function (p) hw(alpha, beta, p)$SSE
                gamma <- optimize(error, lower = 0, upper = 1)$minimum
            }
        }
    } else {
        ## !optimize gamma
        if (is.null(alpha)) {
            ## optimize alpha
            if (is.null(beta)) {
                ## optimize beta
                ## --> optimize alpha and beta
                error <- function (p) hw(p[1L], p[2L], gamma)$SSE
                sol   <- optim(c(optim.start["alpha"], optim.start["beta"]),
                               error, method = "L-BFGS-B",
                               lower = c(0, 0), upper = c(1, 1),
                               control = optim.control)
                if(sol$convergence || any(sol$par < 0 | sol$par > 1)) {
                    if (sol$convergence > 50) {
                        warning(gettextf("optimization difficulties: %s",
                                         sol$message), domain = NA)
                    } else stop("optimization failure")
                }
                alpha <- sol$par[1L]
                beta  <- sol$par[2L]
            } else {
                ## !optimize beta
                ## --> optimize alpha
                error <- function (p) hw(p, beta, gamma)$SSE
                alpha <- optimize(error, lower = 0, upper = 1)$minimum
            }
        } else {
            ## !optimize alpha
            if(is.null(beta)) {
                ## optimize beta
                ## --> optimize beta
                error <- function (p) hw(alpha, p, gamma)$SSE
                beta <- optimize(error, lower = 0, upper = 1)$minimum
            } ## else optimize nothing!
        }
    }

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
