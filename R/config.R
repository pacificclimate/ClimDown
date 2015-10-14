.onLoad <- function(libname, pkgname) {
    op <- options()
    cd.options <- list(
        max.GB=1,
        trimmed.mean=0,
        delta.days=45,
        n.analogues=30,
        obs.ca.years=1951:2005,
        tol=0.1,
        expon=0.5,
        mc.cores=4
    )

    toset <- !(names(cd.options) %in% names(op))

    if(any(toset)) options(cd.options[toset])
}

target.units <- c(tasmax='celsius', tasmin='celsius', pr='mm day-1')

# Takes a vector length and chunk size
# returns a list of (start, stop, length)
chunk.indices <- function(total.size, chunk.size) {
  lapply(
    split(1:total.size, ceiling(1:total.size / chunk.size)),
    function(x) {c('start'=min(x), 'stop'=max(x), 'length'=length(x))}
    )
}

optimal.chunk.size <- function(n.elements, max.GB=getOption('max.GB')) {
  # 8 byte numerics
  floor(max.GB * 2 ** 30 / 8 / n.elements)
}
