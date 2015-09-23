#!/usr/bin/env Rscript

library(ClimDown)

usage <- function() {
    print("Usage: Rscript BCCA.R [gcm_file] [obs_file] [output_file] [variable_name]")
}

# "--args gcm.file='${gcm.file}' obs.file='${obs.file}' output.file='${output.file}' varid='${varid}'"
args <- list(commandArgs(trailingOnly=TRUE))
if (length(args) != 4) {
    usage()
    quit(status=1)
}

names(args) <- c('gcm.file', 'obs.file', 'output.file', 'varid')

gcm <- bias.correct.dqm.netcdf(gcm.file, obs.file, varid)
create.analogues(gcm.file, obs.file, output.file, varid)
