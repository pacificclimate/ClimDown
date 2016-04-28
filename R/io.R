CD_ncvar_get <- function(nc, varid=NA, start=NA, count=NA, verbose=FALSE,
                         signedbyte=TRUE) {
    x <- ncvar_get(nc, varid, start, count, verbose, signedbyte, collapse_degen=FALSE)

    if (var.id == 'pr') {
        x[x < 0] <- 0
    }

    units <- ncatt_get(nc, varid, 'units')$value
    ud.convert(x, units, target.units[varid])
}
