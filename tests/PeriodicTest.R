library(RobustHoltWinters)

test.periodic <- function () {
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
