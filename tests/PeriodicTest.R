library(RobustHoltWinters)

test.periodic <- function () {
    checkAbsLessThanOrEqual <- function(expected, actual, tolerance = 0.001) {
        checkTrue(abs(actual) <= expected + tolerance, paste(actual, 'should be, in absolute value less than', expected))
    }

    synthetic <- function (x) cos(2 * pi * x / 7)

    series <- ts(frequency = 7, synthetic(0:500))
    model <- RobustHoltWinters(series)

    checkEqualsNumeric(0.0, model$coefficients['a'], "Level data should be nearly 0", 0.001)
    checkEqualsNumeric(0.0, model$coefficients['b'], "Trend data should be nearly 0", 0.001)

    checkEqualsNumeric(1.0, model$coefficients['s4'], "Seasonal coefficient should be 1.0", 0.001)

    print(model$coefficients)
    for (i in 1:7) {
        checkAbsLessThanOrEqual(1.0, model$coefficients[paste('s', i, sep='')])
    }
}
