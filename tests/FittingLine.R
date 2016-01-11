library(RobustHoltWinters)

test.constant <- function () {
    series <- ts(frequency = 7, 2 * 0:50)
    model <- RobustHoltWinters(series)

    checkEqualsNumeric(100.0, model$coefficients['a'], "Deseasonalized series should be 0.0, equal to the constant value", 0.001)
    checkEqualsNumeric(2.0, model$coefficients['b'], "Trend be 2.0, equal to the angular coefficient", 0.001)
}
