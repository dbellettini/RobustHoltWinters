library(RobustHoltWinters)

trainingLength <- 100
validationLength <- 1
dataSetLength <- trainingLength + validationLength
iterations <- 100

check <- function (noiseGenerator) {
    allErrors = NULL
    for (i in seq(1, iterations)) {
        tryCatch({
            errors <- simulate.timeseries(noiseGenerator(dataSetLength), c(HoltWinters, RobustHoltWinters))
            allErrors = rbind(allErrors, errors)
        })
    }

    msfes <- NULL
    for (i in seq(1,2)) {
        msfes <- c(msfes, msfe(allErrors[,i]))
    }

    return (msfes)
}

test.cleanData <- function () {
    return(check(function (length) rep(0, length)))
}

test.symmetricOutliers <- function() {
    return(check(function (length) rbinom(length, 1, 0.05) * rnorm(length, 0, 20)))
}

simulate.timeseries <- function (noise, algorithms) {
    data <- arima.sim(list(1,1,1), n = dataSetLength)
    data <- data + noise

    synthetic <- function (x) data[x]
    validation <- synthetic((trainingLength+1):(dataSetLength))

    training <- ts(frequency = 7, synthetic(1:trainingLength))

    errors = NULL
    for (algorithm in algorithms) {
        model <- algorithm(training)
        predictedFutureAmounts <- predict(model, length(validation))
        predicted <- as.numeric(predictedFutureAmounts[,'fit'])

        errors = c(errors, validation[1] - predicted[1])
    }

    # returns errors of 1 step ahead prediction
    return (errors)
}
