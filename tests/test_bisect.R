regrid.one.dim <- ClimDown:::regrid.one.dim

test.regrid.one.dim.monotonic <- function() {
    ## A normal case, with both sequences monotonically increasing
    x <- regrid.one.dim(seq(10, 40, 10), seq(25, 34, 1))
    checkEquals(x, c(2, rep(3, 9)))
}

test.regrid.one.dim.decreasing.fine <- function() {
    ## Fine is decreasing
    x <- regrid.one.dim(seq(10, 100, 10), seq(36, 27, -1))
    checkEquals(x, c(4, rep(3, 9)))
}

test.regrid.one.dim.decreasing.coarse <- function() {
    ## Coarse is decreasing
    x <- regrid.one.dim(seq(100, 10, -10), seq(25, 34, 1))
    checkEquals(x, c(9, rep(8, 9)))
}

test.regrid.one.dim.both.decreasing <- function() {
    ## Both are decreasing
    x <- regrid.one.dim(seq(100, 10, -10), seq(99, 90, -1))
    checkEquals(x, c(rep(1, 4), rep(2, 6)))
}
