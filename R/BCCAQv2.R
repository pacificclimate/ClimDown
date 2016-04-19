#' @title Wrapper function for the entire BCCAQv2 downscaling method
#'
#' @description This function does everything
#'
#' @param gcm.file Filename of GCM simulations in NetCDF format
#' @param obs.file Filename of high-res gridded historical observations
#' @param out.file The file to create (or overwrite) with the final BCCAQv2 NetCDF output
#' @param varname Name of the NetCDF variable to downscale (e.g. 'tasmax')
#'
#' @export
bccaqv2.netcdf.wrapper <- function(gcm.file, obs.file, out.file, varname='tasmax') {
    ptm <- proc.time()

    bcci.file <- tempfile(fileext='.nc')
    bcci.netcdf.wrapper(gcm.file, obs.file, bcci.file, varname)
    qpqm.file <- tempfile(fileext='.nc')
    qpqm.netcdf.wrapper(obs.file, bcci.file, qpqm.file, varname)
    analogues <- bcca.netcdf.wrapper(gcm.file, obs.file, varname)
    qdm.netcdf.wrapper(qpqm.file, obs.file, analogues, out.file, varname)
    unlink(qpqm.file)
    unlink(bcci.file)

    print('Elapsed time')
    print(proc.time() - ptm)
}
