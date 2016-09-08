#' @title Wrapper function for the entire BCCAQ downscaling method
#'
#' @description This function does everything
#'
#' @param gcm.file Filename of GCM simulations in NetCDF format
#' @param obs.file Filename of high-res gridded historical observations
#' @param out.file The file to create (or overwrite) with the final BCCAQ NetCDF output
#' @param varname Name of the NetCDF variable to downscale (e.g. 'tasmax')
#'
#' @references Werner, A. T., & Cannon, A. J. (2016). Hydrologic extremes  - an intercomparison of multiple gridded statistical downscaling methods. Hydrology and Earth System Sciences, 20(4), 1483-1508. doi: 10.5194/hess-20-1483-2016
#' @export
bccaq.netcdf.wrapper <- function(gcm.file, obs.file, out.file, varname='tasmax') {
    ptm <- proc.time()

    ci.file <- tempfile(fileext='.nc')
    ci.netcdf.wrapper(gcm.file, obs.file, ci.file, varname)
    qdm.file <- tempfile(fileext='.nc')
    qdm.netcdf.wrapper(obs.file, ci.file, qdm.file, varname)
    unlink(ci.file)
    analogues <- ca.netcdf.wrapper(gcm.file, obs.file, varname)
    rerank.netcdf.wrapper(qdm.file, obs.file, analogues, out.file, varname)
    unlink(qdm.file)

    print('Elapsed time')
    print(proc.time() - ptm)
}
