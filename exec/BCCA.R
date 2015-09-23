#!/usr/bin/env Rscript

library(ClimDown)

usage <- function() {
    print("Usage: Rscript BCCA.R [gcm_file] [obs_file] [output_file] [variable_name]")
}

args <- list(commandArgs(trailingOnly=TRUE))
if (length(args) != 4) {
    usage()
    quit(status=1)
}

names(args) <- c('gcm.file', 'obs.file', 'output.file', 'varid')
attach(args)

bcca.netcdf.wrapper(gcm.file, obs.file, output.file, varid)
