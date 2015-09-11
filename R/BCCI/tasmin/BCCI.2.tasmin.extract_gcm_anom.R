##******************************************************************************
# Bias Corrected Climate Imprint (BCCI) downscaling algorithm
# Alex Cannon (acannon@uvic.ca)
# *All wind related lines of code - commented out by Arelia Werner (wernera@uvic.ca)
##******************************************************************************
# Read GCM grid and calculate daily anomalies with respect to the GCM monthly
# climatologies
##******************************************************************************
##Modified to only analyze minimum temperature

ptm <- proc.time()

library(ncdf4)

code.dir <- '/home/ssobie/stat.downscaling/code/QPQM/BCCI/'

source(paste(code.dir,'netcdf.calendar.R',sep=''))

args <- commandArgs(trailingOnly=TRUE)
for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
}

config <- paste(code.dir,'BCCI.set.config',sep='')
print(readLines(config))
source(config)

  
  tasmin.nc <- nc_open(tasmin.nc.file)
  
  ##******************************************************************************
  ## Read GCM data
  
  tasmin <- ncvar_get(tasmin.nc, tasmin.var)-273.15

  raw.time <- ncvar_get(tasmin.nc,'time')
  tasmin.dates <- netcdf.calendar(tasmin.nc)
  tasmin.clim.indices <- tasmin.dates[,1] %in% clim.years
  dates <- tasmin.dates

  save(raw.time, file=paste(output.dir, 'tasmin.raw.time', output.suffix, '.RData',
                sep=''))

  save(dates, file=paste(output.dir, 'tasmin.dates.gcm', output.suffix, '.RData',
                sep=''))


  dates <- tasmin.dates
  nc_close(tasmin.nc)
  
################################################################################
# Calculate monthly GCM climatologies

  tasmin.clim <- list()

  for(mn in 1:12){
    tn.mn.clim.indices <- tasmin.clim.indices & (tasmin.dates[,2]==mn)
    tasmin.clim[[mn]] <- apply(tasmin[,,tn.mn.clim.indices], c(1,2), mean)
}

################################################################################
# Calculate daily anomalies with respect to monthly GCM climatologies

  for(i in seq_along(tasmin.dates[,1])){
    cat(tasmin.dates[i,], '\n')
    mn <- tasmin.dates[i,2]
    tasmin[,,i] <- tasmin[,,i]-tasmin.clim[[mn]]
  }
  save(tasmin, file=paste(output.dir, 'tasmin.gcm', output.suffix, '.RData',
                 sep=''))
  rm(tasmin)
  rm(tasmin.dates)
  rm(tasmin.clim)

################################################################################
  ##Using tasmax dates for now



print('Elapsed time')
print(proc.time() - ptm)

rm(list=ls())
##******************************************************************************
##}
