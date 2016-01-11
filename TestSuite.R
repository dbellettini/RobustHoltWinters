library(RUnit)

testSuite <- defineTestSuite("Holt-Winters-Gelper",
                             "tests",
                             ".+")

testResult <- runTestSuite(testSuite)
printTextProtocol(testResult, showDetails = TRUE)
