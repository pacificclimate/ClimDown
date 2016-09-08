#!/usr/bin/env Rscript

library(ClimDown)

usage <- function() {
    print("Usage: Rscript CI.R [gcm_file] [obs_file] [output_file] [variable_name]")
}

args <- as.list(commandArgs(trailingOnly=TRUE))
if (length(args) != 4) {
    usage()
    quit(status=1)
}

names(args) <- c('gcm.file', 'obs.file', 'output.file', 'varid')
attach(args)

ci.netcdf.wrapper(gcm.file, obs.file, output.file, varid)
