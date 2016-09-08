#!/usr/bin/env Rscript

library(ClimDown)

args <- commandArgs(TRUE)

## Help section
usage <- function() {
    cat('
      Bias Correction Constructed Analogues Quantile

      Usage:
      ./QCCAQ.R [gcm_file] [obs_file] [output_file] [variable_name]

      Arguments:
      gcm_file - Filename of GCM simulations in NetCDF format
      obs_file - Filename of high-res gridded historical observations
      output_file - The file to create (or overwrite) with the final BCCAQ NetCDF output
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

bccaq.netcdf.wrapper(gcm.file, obs.file, output.file, varname)
