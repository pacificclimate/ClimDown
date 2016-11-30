#' ClimDown: PCIC's daily Climate Downscaling library
#'
#' This package includes PCIC's routines and techniques for
#' downscaling coarse scale Global Climate Models (GCMs) to fine scale
#' spatial resolution.
#'
#' At present, the package only exports high-level wrapper function
#' that perform each of three downscaling steps: CI, CA, and QDM, as
#' well as one wrapper that runs the entire pipeline: BCCAQ. In
#' general, each wrapper simply takes four arguments: GCM file,
#' gridded observation file, output file, and variable name. However,
#' see the specific function documentation for specifics.
#'
#' The package also provides five wrapper scripts that allow the user
#' to run each step from the command line (plus the whole pipeline)
#' with Rscript.
#'
#' @name ClimDown
#' @aliases ClimDown-package
#' @docType package
#' @references Werner, A. T., & Cannon, A. J. (2016). Hydrologic extremes  - an intercomparison of multiple gridded statistical downscaling methods. Hydrology and Earth System Sciences, 20(4), 1483-1508. doi: 10.5194/hess-20-1483-2016
#' @keywords climate downscaling
#' @import PCICt udunits2 ncdf4 fields foreach seas abind stats
NULL

#' @title User-configurable options
#'
#' @description \pkg{ClimDown} has a number of global options with sensible defaults
#' pre-set. Options must be configured using R's \code{\link[base]{options}}
#' function. They can be set by using an .Rprofile file, or on the R
#' session prompt prior to executing any \pkg{ClimDown} wrapper functions.
#'
#' @param max.GB An \emph{approximate} measure of how much RAM to use in
#'     the chunk I/O loop. In reality, R does a lot of copying data,
#'     so if you have a firm RAM threshold, it's best to set this to
#'     about 1/3 to 1/4 of what you want the high-water mark to
#'     be. (default=1)
#' @param trimmed.mean Undocumented and not recommended to
#'     change. (default=0)
#' @param delta.days Undocumented and not recommended to
#'     change. (default=45)
#' @param n.analogues The number of temporal analogues that the CA
#'     algorithm will search for and match. The higher this number,
#'     the longer the execution time of the reordering
#'     step. (default=30)
#' @param calibration.start A POSIXct object that defines the
#'     beginning of the calibration period. \cr
#'     (default=\code{as.POSIXct('1971-01-01', tz='GMT')})
#' @param calibration.end A POSIXct object that defines the end of the
#'     calibration period. \cr
#'     (default=\code{as.POSIXct('2005-12-31', tz='GMT')})
#' @param tol Undocumented and not recommended to
#'     change. (default=0.1)
#' @param expon Undocumented and not recommended to
#'     change. (default=0.5)
#'
#' @param multiyear Undocumented and not recommended to
#'     change. (default=TRUE)
#' @param expand.multiyear Undocumented and not recommended to
#'     change. (default=TRUE)
#' @param multiyear.window.length Undocumented and not recommended to
#'     change. (default=30)
#' @param trace Undocumented and not recommended to
#'     change. (default=0.005)
#' @param jitter.factor Undocumented and not recommended to
#'     change. (default=0.01)
#' @param tau Undocumented and not recommended to change. \cr
#'     (default=\code{list(pr=1001, tasmax=101, tasmin=101)})
#' @param seasonal Undocumented and not recommended to change. \cr
#'     (default=\code{list(pr=TRUE, tasmax=FALSE, tasmin=FALSE)})
#' @param ratio Undocumented and not recommended to change. \cr
#'     (default=\code{list(pr=TRUE, tasmax=FALSE, tasmin=FALSE)})
#'
#' @param check.units A boolean value that determines whether to check
#'     the input units and convert them to the target output
#'     units. The \emph{safe} option is to leave this set to TRUE, but
#'     if know for sure that the units match, modest performance gains
#'     can be made by not checking and performing this
#'     conversion. (default=TRUE)
#'
#' @param check.neg.precip A boolean value that determines whether to
#'     check for and eliminate negative precipitation values. Like
#'     check.units, modest performance gains can be achieved if you
#'     are certain that no negative precipitation values
#'     exist. (default=TRUE)
#'
#' @param target.units A list containing the units that should be used
#'     for the output file \cr
#'     (default=\code{c(tasmax='celsius', tasmin='celsius', pr='kg m-2 d-1')})
#'
#' @name options
NULL

#' @title Parallelization
#'
#' @description \pkg{ClimDown} uses the \pkg{foreach} package to
#' support parallelization where possible. In general, most of
#' \pkg{ClimDown}'s computations contain an outer I/O chunk loop and
#' an inner computational loop. The outer chunk loop serially reads in
#' as much data from input files as it feasibly can, and then the
#' inner loop will execute in parallel if the user has configured a
#' parallel engine.
#'
#' Users can configure a parallel engine before execution using either
#' the \pkg{doParallel} or \pkg{doMPI} packages.
#'
#' If the user does not configure a parallel backend, the inner loops
#' will simply run serially and issue a warning.
#'
#' @examples
#' \dontrun{
#' library(doParallel)
#' registerDoParallel(cores=4)
#' bccaq.netcdf.wrapper(...)
#' stopImplicitCluster()}
#'
#' @seealso \pkg{doParallel} and \pkg{doMPI}
#' @name parallelization
NULL
