#!/usr/bin/env Rscript

args <- commandArgs(TRUE)

## Help section
usage <- function() {
    cat('
      Bias Correction Constructed Analogues Quantile version 2

      Usage:
      ./QCCAQv2.R [gcm_file] [obs_file] [output_file] [variable_name]

      Arguments:
      gcm_file - Filename of GCM simulations in NetCDF format
      obs_file - Filename of high-res gridded historical observations
      output_file - The file to create (or overwrite) with the final BCCAQv2 NetCDF output
      varname - The name of the variable to downscale'
    )

    q(save='no')
}

args <- as.list(commandArgs(trailingOnly=TRUE))
if (length(args) != 4) {
    usage()
    quit(status=1)
}

names(args) <- c('gcm.file', 'obs.file', 'output.file', 'varname')
attach(args)

bccaqv2.netcdf.wrapper(gcm.file, obs.file, output.file, varname)
