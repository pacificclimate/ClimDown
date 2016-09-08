#!/usr/bin/env Rscript

library(ClimDown)

args <- commandArgs(TRUE)

## Help section
usage <- function() {
    cat('
      Quantile Reranking

      Usage:
      ./rerank.R [qdm_file] [obs_file] [analogues_file] [output_file] [variable_name]

      Arguments:
      qdm_file - The output file from the QDM script
      obs_file - Filename of high-res gridded historical observations
      analogues_file - The output file from the BCCA script
      output_file - The file to create (or overwrite) with the ...
      varname - The name of the variable to downscale

      All files (save for the analogues_file) should have the same spatial domain.\n\n'
    )

    q(save='no')
}

args <- as.list(commandArgs(trailingOnly=TRUE))
if (length(args) != 5) {
    usage()
    quit(status=1)
}

names(args) <- c('qdm.file', 'obs.file', 'analogues.file', 'output.file', 'varid')
attach(args)

load(analogues.file)
rerank.netcdf.wrapper(gcm.file, bcci.file, obs.file, analogues, output.file, varid)
