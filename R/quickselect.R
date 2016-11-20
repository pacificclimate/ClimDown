introselect.k.lowest <- function(thelist, k) {
    kth.lowest <- introselect(thelist, k)
    which(thelist <= kth.lowest)[1:k]
}
