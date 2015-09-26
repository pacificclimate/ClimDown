# Implementation of the quickselect algorithm for selecting
# the kth lowest value in an unsorted list
# James Hiebert <hiebert@uvic.ca>

# Note the swapping makes for kind of awkward/messy code, but one really
# needs to avoid doing any dynamic memory allocation in this core code
partition <- function(thelist, left, right, pivot.i) {
    val <- thelist[pivot.i]
    x <- thelist[right]
    thelist[right] <- thelist[pivot.i]
    thelist[pivot.i] <- x

    storeIndex <- left
    for (i in left:(right-1)) {
        if (thelist[i] < val) {
            x <- thelist[storeIndex]
            thelist[storeIndex] <- thelist[i]
            thelist[i] <- x
            storeIndex <- storeIndex + 1
        }
    }
    x <- thelist[right]
    thelist[right] <- thelist[storeIndex]
    thelist[storeIndex] <- x

    list(list=thelist, storeIndex=storeIndex)
}

# Return the value of the kth lowest element in thelist
# in between elements indexed by left and right
quickselect <- function(thelist, k, left, right) {
    # Base case: return the single element if it's the only one
    if (left == right) {
        return(thelist[left])
    }

    # Choose a pivot point at random
    pivot.i <- sample(left:right, 1)

    rv <- partition(thelist, left, right, pivot.i)
    thelist <- rv$list
    pivot.i <- rv$storeIndex
    # The pivot is in its final sorted position
    if (pivot.i == k) {
        thelist[k]
    }
    # Our search value is in the left hand branch
    else if (k < pivot.i) {
        quickselect(thelist, k, left, pivot.i - 1)
    }
    # Our search value is in the right hand branch
    else {
        quickselect(thelist, k, pivot.i + 1, right)
    }
}

# Return the indices of all k elements that are the lowest elements in thelist
quickselect.k.lowest <- function(thelist, k, left=1, right=NULL) {
    if (is.null(right)) {
        right <- length(thelist)
    }
    kth.lowest <- quickselect(thelist, k, left, right)
    # In the unlikely possibility of ties, simply take the first k elements
    which(thelist <= kth.lowest)[1:k]
}

# In R this is the simplest way to select the kth lowest value in a list, but
# rank() has to actually sort all values of the list which means it can be no better than
# O(n log n). The advantage, however is that it *does* have a low constant cost
# so it turns out to be faster for smaller values of n
slowselect.k.lowest <- function(thelist, k) {
    which(rank(thelist, ties.method='random') <= k)
}

# Decide which algorithm to use based on the size of the input
# One million seems to be the cutoff where the quickselect algorithm is significantly
# faster than the naive ordering algorithm
quickest.select <- function(l, k) {
    if (length(l) < 1000000) {
        slowselect.k.lowest(l, k)
    } else {
        quickselect.k.lowest(l, k)
    }       
}
