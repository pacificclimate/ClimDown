#!/usr/bin/env Rscript

args <- commandArgs(TRUE)

## Help section
usage <- function() {
    cat('
      Quantile Delta Mapping???

      Usage:
      ./QDM.R [qpqm_file] [bcci_file] [analogues_file] [output_file] [variable_name]

      Arguments:
      qpqm_file - The output file from the QPQM script
      bcci_file - The output file from the BCCI script
      obs_file - Filename of high-res gridded historical observations
      analogues_file - The output file from the BCCA script
      output_file - The file to create (or overwrite) with the ...
      varname - The name of the variable to downscale

      All files (save for the analogues_file) should have the same spatial domain.\n\n'
    )

    q(save='no')
}

args <- as.list(commandArgs(trailingOnly=TRUE))
if (length(args) != 6) {
    usage()
    quit(status=1)
}

names(args) <- c('qpqm.file', 'bcci.file', 'obs.file', 'analogues.file', 'output.file', 'varid')
attach(args)

load(analogues.file)
qdm.netcdf.wrapper(gcm.file, bcci.file, obs.file, analogues, output.file, varid)
