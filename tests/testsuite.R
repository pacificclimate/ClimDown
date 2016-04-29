library(RUnit)
library(ncdf4)

if (require("RUnit", quietly=TRUE)) {
    wd <- getwd()
    loadNamespace("ClimDown") ## Ensure default options are set
    testsuite <- defineTestSuite("ClimDown", dirs=wd, testFileRegexp = "^test_.+.R$", testFuncRegexp = "^test.+")
    test.result <- runTestSuite(testsuite, useOwnErrorHandler=F)
    printTextProtocol(test.result)
    stopifnot(test.result$ClimDown$nFail == 0 && test.result$ClimDown$nErr == 0)
}
