##******************************************************************************
# Bias Corrected Climate Imprint (BCCI) downscaling algorithm
# Alex Cannon (acannon@uvic.ca)
# *All wind related lines of code - commented out by Arelia Werner (wernera@uvic.ca)
##******************************************************************************
# Read fine-scale grid and calculate monthly climatologies
##******************************************************************************

library(ncdf4)
library(RNetCDF)
source('netcdf.calendar.R')


config <- commandArgs(trailingOnly=TRUE)
##
if(length(config)==0) config <- 'BCCI.set.config'
print(readLines(config))
source(config)

##  source(config) ##Should be a unique configure file
  nc <- nc_open(nc.obs.file)

##******************************************************************************

  obs.lon <- ncvar_get(nc, 'lon')
  obs.lat <- ncvar_get(nc, 'lat')
  lonlat.clim <- list(lon=obs.lon, lat=obs.lat)
  obs.time <- netcdf.calendar(nc)

################################################################################
  ## Calculate monthly fine-scale climatologies

  n.lon <- length(obs.lon)
  n.lat <- length(obs.lat)
  n.time <- nrow(obs.time)
  
  clim.indices <- obs.time[,1] %in% clim.years
  st <- head(which(clim.indices),1) - 1
  en <- tail(which(clim.indices),1) - 1

##Only needed if the subset and averaging has not been completed:
system(paste('ncks -O -d time,',st,',',en,' ',nc.obs.file,' ',nc.obs.sub,sep=''))
system(paste('cdo -s ymonmean ',nc.obs.sub,' ',nc.obs.clim,sep=''))

nc.clim <- nc_open(nc.obs.clim)

pr.clim <- tasmax.clim <- tasmin.clim <- list()
  ## <- wind.clim 
  for(mn in 1:12){
    pr.clim[[mn]] <- ncvar_get(nc.clim,'pr',start=c(1,1,mn),count=c(-1,-1,1))
    tasmax.clim[[mn]] <- ncvar_get(nc.clim,'tasmax',start=c(1,1,mn),count=c(-1,-1,1))
    tasmin.clim[[mn]] <- ncvar_get(nc.clim,'tasmax',start=c(1,1,mn),count=c(-1,-1,1))
  }

  




  nc_close(nc)
  nc_close(nc.clim)
################################################################################
  
  save(lonlat.clim, file=paste(output.dir, 'lonlat.obs.clim', output.suffix,
                      '.RData', sep=''))
  save(pr.clim, file=paste(output.dir, 'pr.obs.clim', output.suffix,
                  '.RData', sep=''))
  save(tasmax.clim, file=paste(output.dir, 'tasmax.obs.clim', output.suffix,
                      '.RData', sep=''))
  save(tasmin.clim, file=paste(output.dir, 'tasmin.obs.clim', output.suffix,
                      '.RData', sep=''))
  ## save(wind.clim, file=paste(output.dir, 'wind.obs.clim', output.suffix,
  ## '.RData', sep=''))
  rm(list=ls())
##******************************************************************************
##}
