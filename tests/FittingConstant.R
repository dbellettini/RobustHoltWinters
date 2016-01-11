library(RobustHoltWinters)

test.deseasonalized <- function () {
    series <- ts(frequency = 7, rep(10, 50))
    model <- RobustHoltWinters(series)

    checkEqualsNumeric(10.0, model$coefficients['a'], "Deseasonalized series should be 10.0, equal to the constant value", 0.001)
}
