#' @title Wrapper function for the entire BCCAQ downscaling method
#'
#' @description BCCAQ is a hybrid downscaling method that combines
#' outputs from Climate Analogues (CA) and quantile mapping at the
#' fine-scale resolution.  First, the CA and climate imprint (CI)
#' plus quantile delta mapping (QDM) algorithms are run
#' independently. BCCAQ then combines outputs from the two by
#' taking the daily QDM outputs at each fine-scale grid point and
#' reordering them within a given month according to the daily
#' CA ranks, i.e., using a form of Empirical Copula
#' Coupling.
#'
#' The combination mitigates some potential issues with
#' the separate algorithms. First, because the optimal weights
#' used to combine the analogues in BCCA are derived on a
#' day-by-day basis, without reference to the full historical data
#' set, the algorithm may fail to reproduce long-term trends from
#' the climate model. Second, the CI/QDM bias correction step
#' fixes precipitation "drizzle" and other residual biases caused
#' by the linear combination of daily fields from CA. Third,
#' reordering data for each fine-scale grid point within a month
#' effectively breaks the overly smooth representation of sub
#' grid-scale spatial variability inherited from CI/QDM, thereby
#' resulting in a more accurate representation of event-scale
#' spatial gradients; this also prevents the downscaled outputs
#' from drifting too far from the climate model's long-term trend.
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
