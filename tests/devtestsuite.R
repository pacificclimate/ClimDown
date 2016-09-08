library(RUnit)
library(ncdf4)
library(udunits2)
library(PCICt)
library(fields)
library(foreach)

if (require("RUnit", quietly=TRUE)) {
    wd <- getwd()
    devtools::load_all("..") ## Ensure default options are set
    testsuite <- defineTestSuite("ClimDown_dev", dirs=wd, testFileRegexp = "^test.+.R$", testFuncRegexp = "^test.+")
    test.result <- runTestSuite(testsuite, useOwnErrorHandler=F)
    printTextProtocol(test.result)
    stopifnot(test.result$ClimDown$nFail == 0 && test.result$ClimDown$nErr == 0)
}
