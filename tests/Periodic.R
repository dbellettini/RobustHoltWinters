library(RobustHoltWinters)

test.periodic <- function () {
    synthetic <- function (x) cos(2 * pi * x / 7)

    series <- ts(frequency = 7, synthetic(0:500))
    model <- RobustHoltWinters(series)

    print(model)
    checkEqualsNumeric(1.0, model$coefficients['s4'], "Seasonal coefficient should be 1.0", 0.001)
}
