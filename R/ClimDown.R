#' ClimDown: PCIC's daily Climate Downscaling library
#'
#' This package includes PCIC's routines and techniques for
#' downscaling coarse scale Global Climate Models (GCMs) to fine scale
#' spatial resolution.
#'
#' At present, the package only exports high-level wrapper function
#' that perform each of three downscaling steps: BCCI, BCCA, and
#' QPQM. Each wrapper simply takes four arguments: GCM file, gridded
#' observation file, output file, and variable name.
#'
#' The package also provides three wrapper scripts that allow the user
#' to run each step from the command line with Rscript.
#'
#' @name ClimDown
#' @aliases ClimDown-package
#' @docType package
#' @references Alex J. Cannon, Stephen R. Sobie, and Trevor Q. Murdock, 2015: Bias Correction of GCM Precipitation by Quantile Mapping: How Well Do Methods Preserve Changes in Quantiles and Extremes?. J. Climate, 28, 6938–6959.
#' @keywords climate downscaling
#' @import PCICt udunits2 ncdf4 fields foreach seas
NULL
