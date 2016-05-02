.setUp <- function() {
    bcci.file <<- tempfile(fileext='.nc')
    ClimDown::bcci.netcdf.wrapper('./tiny_gcm.nc', './tiny_obs.nc', bcci.file, 'tasmax')
}

.tearDown <- function() {
    unlink(bcci.file)
}

test.qpqm <- function() {
    out.nc <- tempfile(fileext='.nc')
    options(
        calibration.end=as.POSIXct('1972-12-31', tz='GMT'),
        cend=as.POSIXct('1972-12-31', tz='GMT')
    )
    ClimDown::qpqm.netcdf.wrapper('./tiny_obs.nc', bcci.file, out.nc, 'tasmax')
    unlink(out.nc)
    checkTrue(TRUE)
}
