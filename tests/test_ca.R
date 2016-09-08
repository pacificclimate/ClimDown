test.ca.netcdf.wrapper <- function() {
    options(
        calibration.end=as.POSIXct('1972-12-31', tz='GMT')
    )
    analogues <- ClimDown::ca.netcdf.wrapper('./tiny_gcm.nc', './tiny_obs.nc')
    checkEquals(names(analogues), c("indices", "weights"))
    checkTrue(is.integer(analogues$indices[[1]]))
    checkTrue(is.numeric(analogues$weights[[1]]))
}

test.ca.netcdf.wrapper.bad.times <- function() {
    options(calibration.end=as.POSIXct('2005-12-31', tz='GMT'))
    checkException(ClimDown::ca.netcdf.wrapper('./tiny_gcm.nc', './tiny_obs.nc'))
}
