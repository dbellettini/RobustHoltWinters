source('tests/Setup.R')

library(RUnit)

testSuite <- defineTestSuite("Holt-Winters-Gelper",
                             "tests",
                             ".*Outlier.+")

testResult <- runTestSuite(testSuite)
printTextProtocol(testResult, showDetails = TRUE)
