.setUp <- function() {
    options(
        calibration.end=as.POSIXct('1972-12-31', tz='GMT')
    )
    bcci.file <- tempfile(fileext='.nc')
    ClimDown::bcci.netcdf.wrapper('./tiny_gcm.nc', './tiny_obs.nc', bcci.file, 'tasmax')
    qpqm.file <<- tempfile(fileext='.nc')
    ClimDown::qpqm.netcdf.wrapper('./tiny_obs.nc', bcci.file, qpqm.file, 'tasmax')
    unlink(bcci.file)
    analogues <<- ClimDown::bcca.netcdf.wrapper('./tiny_gcm.nc', './tiny_obs.nc')
}

.tearDown <- function() {
    unlink(qpqm.file)
}

test.qdm <- function() {
    out.file <- tempfile(fileext='.nc')
    ClimDown::qdm.netcdf.wrapper(qpqm.file, './tiny_obs.nc', analogues, out.file, varname='tasmax')
    unlink(out.file)
    checkTrue(TRUE)
}
