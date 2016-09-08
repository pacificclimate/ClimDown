#!/usr/bin/env Rscript

library(ClimDown)

args <- commandArgs(TRUE)

## Help section
usage <- function() {
    cat('
      Quantile Delta Mapping

      Usage:
      ./QDM.R [obs_file] [gcm_file] [output_file] [variable_name]

      Arguments:
      obs_file - Gridded historical observations
      gcm_file - GCM simulations interpolated to the obs_file grid
      output_file - The file to create (or overwrite) with the bias corrected outputs
      varname - The name of the variable to downscale

      All files should have the same spatial domain.\n\n'
        )
    
    q(save='no')
}

args <- as.list(commandArgs(trailingOnly=TRUE))
if (length(args) != 4) {
    usage()
    quit(status=1)
}

names(args) <- c('gcm.file', 'obs.file', 'output.file', 'varid')
attach(args)

qdm.netcdf.wrapper(gcm.file, obs.file, output.file, varid)
