test.introselect <- function() {
    x <- c(2, 4, 1, 5, 3, 6, 9, 8, 7, 10)
    for (i in 1:10) {
        rv <- introselect(x, i)
        checkEquals(rv, i)
    }
}
