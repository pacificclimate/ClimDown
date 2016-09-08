.onLoad <- function(libname, pkgname) {
    op <- options()
    cd.options <- list(
        # Computation options
        max.GB=1,
        # CA options
        trimmed.mean=0,
        delta.days=45,
        n.analogues=30,
        calibration.start=as.POSIXct('1971-01-01', tz='GMT'),
        calibration.end=as.POSIXct('2005-12-31', tz='GMT'),
        tol=0.1,
        expon=0.5,
        # QDM options
        multiyear=TRUE,
        expand.multiyear=TRUE,
        multiyear.window.length=30,
        trace=0.005,
        jitter.factor=0.01,
        tau=list(pr=1001, tasmax=101, tasmin=101),
        seasonal=list(pr=TRUE, tasmax=FALSE, tasmin=FALSE),
        ratio=list(pr=TRUE, tasmax=FALSE, tasmin=FALSE),
        # Data processing options
        check.units=TRUE,
        check.neg.precip=TRUE,
        target.units=c(tasmax='celsius', tasmin='celsius', pr='kg m-2 d-1') ##pr='mm day-1')
    )

    toset <- !(names(cd.options) %in% names(op))

    if(any(toset)) options(cd.options[toset])
}


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

# Takes a vector of PCICt dates and chunk.size and splits the vector into chunks
# that are *approximately* of that size, but only break on the boundaries
# between months
chunk.month.factor <- function(t, chunk.size) {
    time.factor <- factor(format(t, '%Y-%m'))
    chunk.factor <- factor(ceiling(1:length(t) / chunk.size))
    f <- interaction(time.factor, chunk.factor, drop=T)

    # Do two passes across the levels of the factor:
    # The first pass merges a month across chunk boundaries
    nl <- nlevels(f)
    new.levels <- mapply(
        function(prev, this) {
            prev.month <- strsplit(prev, '.', fixed=T)[[1]][1]
            this.month <- strsplit(this, '.', fixed=T)[[1]][1]
            if (prev.month == this.month) {
                prev
            }
            else {
                this
            }
        },
        levels(f)[1:nl-1], levels(f)[2:nl]
        )
    # The second pass eliminates the months from the factor
    levels(f) <- gsub('.*\\.(.*)', '\\1', c(levels(f)[1], new.levels))
    f
}
