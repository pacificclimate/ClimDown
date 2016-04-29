test.bcca.netcdf.wrapper <- function() {
    analogues <- ClimDown::bcca.netcdf.wrapper('./tiny_gcm.nc', './tiny_obs.nc')
    checkEquals(names(analogues), c("indices", "weights"))
    checkTrue(is.integer(analogues$indices[[1]]))
    checkTrue(is.numeric(analogues$weights[[1]]))
}
