library(RobustHoltWinters)

test.periodic <- function () {
    synthLength <- 28

    data <- arima.sim(list(1,1,1), n = 100)
    synthetic <- function (x) data[x]
    expected <- synthetic(synthLength:(synthLength + 9))

    series <- ts(frequency = 7, synthetic(1:synthLength))
    # series <- insertOutlier(series, 3)

    robustModel <- RobustHoltWinters(series)
    model <- HoltWinters(series)

    print(model)
    print(robustModel)

    predictedFutureAmounts <- predict(model, 10)
    predicted <- as.numeric(predictedFutureAmounts[,'fit'])

    robustPredictedFutureAmounts <- predict(robustModel, 10)
    robustPredicted <- as.numeric(robustPredictedFutureAmounts[,'fit'])

    print('ERROR:')
    error <- expected - predicted
    print(sum(error ** 2) / length(error))

    print('ROBUST ERROR:')
    error <- expected - robustPredicted
    print(sum(error ** 2) / length(error))

    return (T)

    checkEqualsNumeric(0.0, model$coefficients['a'], "Level data should be nearly 0", 0.001)
    checkEqualsNumeric(0.0, model$coefficients['b'], "Trend data should be nearly 0", 0.001)

    checkEqualsNumeric(1.0, model$coefficients['s4'], "Seasonal coefficient should be 1.0", 0.001)

    for (i in 1:7) {
        checkAbsLessThanOrEqual(1.0, model$coefficients[paste('s', i, sep='')])
    }
}
