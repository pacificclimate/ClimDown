## Perform a binary search on the *sorted* (ascending) vector v
## Return the array index of the element closest to x
find.nearest <- function(x, v) {
    stopifnot(is.monotonic(v))
    if (length(v) == 1) {
        return(1)
    }
    if (length(v) == 2) {
        return(which.min(abs(v - x)))
    }
    mid <- ceiling(length(v) / 2)
    if (x == v[mid]) {
        return(mid)
    } else if (x < v[mid]) {
        return(find.nearest(x, v[1:mid]))
    }
    else {
        return((mid - 1) + find.nearest(x, v[mid:length(v)]))
    }
}

regrid.one.dim <- function(coarse.points, fine.points) {
    if (is.monotonic(coarse.points)) {
        sapply(fine.points, find.nearest, coarse.points)
    } else {
        ## If grid points aren't in sorted order, first save a map of the
        ## sorted array indexes...
        sorted.coarse <- sort.int(coarse.points, index.return=TRUE)
        i <- sapply(fine.points, find.nearest, sorted.coarse$x)
        ## ... then reapply the index map
        sorted.coarse$ix[i]
    }
}

## Take a fine scale (e.g. ANUSPLINE) grid of latitudes and longitudes
## and find the indicies that correspond to a coarse scale (e.g. a GCM) grid
## Since the search is essentially a minimizing distance in 2 dimensions
## We can actually search independently in each dimensions separately (which
## is a huge optimization, making the run time x + y instead of x * y) and
## then reconstruct the indices to create a full grid
regrid.coarse.to.fine <- function(coarse.lats, coarse.lons, fine.lats, fine.lons) {
    xi <- regrid.one.dim(coarse.lons, fine.lons)
    yi <- regrid.one.dim(coarse.lats, fine.lats)

    ## Two dimensional grid of indices
    xi <- matrix(xi, ncol=length(fine.lats), nrow=length(fine.lons), byrow=F)
    yi <- matrix(yi, ncol=length(fine.lats), nrow=length(fine.lons), byrow=T)
    return(list(xi=xi, yi=yi))
}

is.monotonic <- function(x) {
    all(x == cummax(x))
}
