# test suite setup (common functions)
install.packages('.', repos=NULL, type="source")
library(RobustHoltWinters)

checkAbsLessThanOrEqual <- function(expected, actual, tolerance = 0.001) {
    checkTrue(abs(actual) <= expected + tolerance, paste(actual, 'should be, in absolute value less than', expected))
}

insertOutlier <- function(x, moltiplicativeNoise, daysback = 7) {
    i <- length(x) - daysback;
    x[i] <- x[i] * moltiplicativeNoise
    return (x)
}

msfe <- function (error) {
    return (sum(error ** 2) / length(error))
}

performanceEvaluation <- function(model, validation) {
    # TODO: implement performanceEvaluation
    predictedFutureAmounts <- predict(model, length(validation))
    predicted <- as.numeric(predictedFutureAmounts[,'fit'])

    return (msfe(validation, predicted))
}
