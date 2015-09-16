## Perform a binary search on the *sorted* vector v
## Return the array index of the element closest to x
find.nearest <- function(x, v) {
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
    return(sapply(fine.points, find.nearest, coarse.points))
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
