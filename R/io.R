positive_pr <- function(x, varid) {
    if (varid == 'pr') {
        x[x < 0] <- 0
    }
    x
}

CD_ncvar_get <- function(nc, varid=NA, start=NA, count=NA, verbose=FALSE,
                         signedbyte=TRUE, source.units=NULL) {
    x <- ncvar_get(nc, varid, start, count, verbose, signedbyte, collapse_degen=FALSE)
    if (getOption('check.neg.precip')) {
        x <- positive_pr(x, varid)
    }
    if (getOption('check.units')) {
        if (is.null(source.units)) {
            source.units <- ncatt_get(nc, varid, 'units')$value
        }
        ud.convert(x, source.units, getOption('target.units')[varid])
    } else {
        x
    }
}
