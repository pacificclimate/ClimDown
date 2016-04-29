test.bcca.netcdf.wrapper <- function() {
    options(
        calibration.end=as.POSIXct('1972-12-31', tz='GMT')
    )
    analogues <- ClimDown::bcca.netcdf.wrapper('./tiny_gcm.nc', './tiny_obs.nc')
    checkEquals(names(analogues), c("indices", "weights"))
    checkTrue(is.integer(analogues$indices[[1]]))
    checkTrue(is.numeric(analogues$weights[[1]]))
}

test.bcca.netcdf.wrapper.bad.times <- function() {
    options(calibration.end=as.POSIXct('2005-12-31', tz='GMT'))
    checkException(ClimDown::bcca.netcdf.wrapper('./tiny_gcm.nc', './tiny_obs.nc'))
}
