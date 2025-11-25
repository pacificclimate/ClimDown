.setUp <- function() {
    ci.file <<- tempfile(fileext='.nc')
    ClimDown::ci.netcdf.wrapper('./tiny_gcm.nc', './tiny_obs.nc', ci.file)
}

.tearDown <- function() {
    unlink(ci.file)
}

test.qdm <- function() {
    out.nc <- tempfile(fileext='.nc')
    options(
        gcm.varname       = "tasmax",
        obs.varname       = "tasmax",
        calibration.end=as.POSIXct('1972-12-31', tz='GMT'),
        cend=as.POSIXct('1972-12-31', tz='GMT')
    )
    ClimDown::qdm.netcdf.wrapper('./tiny_obs.nc', ci.file, out.nc, 'tasmax')
    unlink(out.nc)
    checkTrue(TRUE)
}
