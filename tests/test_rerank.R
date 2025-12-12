.setUp <- function() {
    options(
        gcm.varname="tasmax",
        obs.varname="tasmax",
        calibration.end=as.POSIXct('1972-12-31', tz='GMT')
    )
    ci.file <- tempfile(fileext='.nc')
    ClimDown::ci.netcdf.wrapper('./tiny_gcm.nc', './tiny_obs.nc', ci.file)
    qdm.file <<- tempfile(fileext='.nc')
    ClimDown::qdm.netcdf.wrapper('./tiny_obs.nc', ci.file, qdm.file, 'tasmax')
    unlink(ci.file)
    analogues <<- ClimDown::ca.netcdf.wrapper('./tiny_gcm.nc', './tiny_obs.nc')
}

.tearDown <- function() {
    unlink(qdm.file)
}

test.rerank <- function() {
    out.file <- tempfile(fileext='.nc')
    ClimDown::rerank.netcdf.wrapper(qdm.file, './tiny_obs.nc', analogues, out.file, varname='tasmax')
    unlink(out.file)
    checkTrue(TRUE)
}
