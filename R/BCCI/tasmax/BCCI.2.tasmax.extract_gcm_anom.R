##******************************************************************************
# Bias Corrected Climate Imprint (BCCI) downscaling algorithm
# Alex Cannon (acannon@uvic.ca)
# *All wind related lines of code - commented out by Arelia Werner (wernera@uvic.ca)
##******************************************************************************
# Read GCM grid and calculate daily anomalies with respect to the GCM monthly
# climatologies
##******************************************************************************
##Modified to only analyze maximum temperature

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

tasmax.nc <- nc_open(tasmax.nc.file)
  
  ##******************************************************************************
  ## Read GCM data
  
  tasmax <- ncvar_get(tasmax.nc, tasmax.var)-273.15
##  lat <- ncvar_get(tasmax.nc, 'lat')
##  lon <- ncvar_get(tasmax.nc, 'lon')
##  lonlat <- list(lon=lon, lat=lat)
##  save(lonlat, file=paste(output.dir, 'lonlat.gcm', output.suffix, '.RData',
##                 sep=''))
##  rm(lonlat)

  raw.time <- ncvar_get(tasmax.nc,'time')
  tasmax.dates <- netcdf.calendar(tasmax.nc)
  tasmax.clim.indices <- tasmax.dates[,1] %in% clim.years
  dates <- tasmax.dates

  save(raw.time, file=paste(output.dir, 'tasmax.raw.time', output.suffix, '.RData',
                sep='')) 

  save(dates, file=paste(output.dir, 'tasmax.dates.gcm', output.suffix, '.RData',
                sep='')) 

  nc_close(tasmax.nc)
##  rm(raw.time)
  
################################################################################
# Calculate monthly GCM climatologies

  tasmax.clim <- list()

  for(mn in 1:12){
    tx.mn.clim.indices <- tasmax.clim.indices & (tasmax.dates[,2]==mn)
    tasmax.clim[[mn]] <- apply(tasmax[,,tx.mn.clim.indices], c(1,2), mean)
  }
################################################################################
# Calculate daily anomalies with respect to monthly GCM climatologies

  for(i in seq_along(tasmax.dates[,1])){
    cat(tasmax.dates[i,], '\n')
    mn <- tasmax.dates[i,2]
    tasmax[,,i] <- tasmax[,,i]-tasmax.clim[[mn]]
  }
  save(tasmax, file=paste(output.dir, 'tasmax.gcm', output.suffix, '.RData',
                 sep=''))
  rm(tasmax)
  rm(tasmax.dates)
  rm(tasmax.clim)	

################################################################################
  ##Using tasmax dates for now

print('Elapsed time')
print(proc.time() - ptm)

#mail("ssobie@uvic.ca","BCCI.2 finished",paste("The BCCI.2 script has finished with ",gcm,"-",rcm,sep=""))

rm(list=ls())
##******************************************************************************
##}
