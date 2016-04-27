library(RUnit)
library(ncdf4)
#devtools::load_all('~/code/git/ClimDown')

if (require("RUnit", quietly=TRUE)) {
    wd <- getwd()
    testsuite <- defineTestSuite("ClimDown", dirs=wd, testFileRegexp = "^test_.+.R$", testFuncRegexp = "^test.+")
    test.result <- runTestSuite(testsuite, useOwnErrorHandler=F)
    printTextProtocol(test.result)
}
