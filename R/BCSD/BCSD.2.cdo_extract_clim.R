##******************************************************************************
# Bias Corrected Spatial Disaggregation (BCSD) downscaling algorithm
# Alex Cannon (acannon@uvic.ca)
##******************************************************************************
# Read fine-scale grids and calculate monthly climatologies
# Calculate min/max values for possible range checking during temporal
# disaggregation.
##******************************************************************************

rm(list=ls())
library(ncdf4)
library(RNetCDF)
source('netcdf.calendar.R')

config <- commandArgs(trailingOnly=TRUE)
if(length(config)==0) config <- 'BCSD.config'
print(readLines(config))
source(config)

nc.gcm <- nc_open(pr.nc.file)
gcm.time <- netcdf.calendar(nc.gcm)
nc_close(nc.gcm)
save(gcm.time, file=paste(output.dir, 'gcm.time', output.suffix,
     '.RData', sep=''))

##******************************************************************************
# Calculate climatological statistics of fine-scale data using cdo

system(paste('cdo -ymonmean -seldate,', min(clim.years), '-01-01T00:00,',
             max(clim.years), '-12-31T23:59 ', nc.obs.file, ' ', output.dir,
             'clim', output.suffix, '.RData', sep=''))
system(paste('cdo -ymonmax -seldate,', min(clim.years), '-01-01T00:00,',
             max(clim.years), '-12-31T23:59 ', nc.obs.file, ' ', output.dir,
             'max', output.suffix, '.RData', sep=''))
system(paste('cdo -ymonmin -seldate,', min(clim.years), '-01-01T00:00,',
             max(clim.years), '-12-31T23:59 ', nc.obs.file, ' ', output.dir,
             'min', output.suffix, '.RData', sep=''))

nc.clim <- nc_open(paste(output.dir, 'clim', output.suffix, '.RData', sep=''))
nc.max <- nc_open(paste(output.dir, 'max', output.suffix, '.RData', sep=''))
nc.min <- nc_open(paste(output.dir, 'min', output.suffix, '.RData', sep=''))

pr.clim.mat <- ncvar_get(nc.clim, 'pr')
tasmax.clim.mat <- ncvar_get(nc.clim, 'tasmax')
tasmin.clim.mat <- ncvar_get(nc.clim, 'tasmin')
wind.clim.mat <- ncvar_get(nc.clim, 'wind')

pr.max.mat <- ncvar_get(nc.max, 'pr')
tasmax.max.mat <- ncvar_get(nc.max, 'tasmax')
tasmin.max.mat <- ncvar_get(nc.max, 'tasmin')
wind.max.mat <- ncvar_get(nc.max, 'wind')

pr.min.mat <- ncvar_get(nc.min, 'pr')
tasmax.min.mat <- ncvar_get(nc.min, 'tasmax')
tasmin.min.mat <- ncvar_get(nc.min, 'tasmin')
wind.min.mat <- ncvar_get(nc.min, 'wind')

nc_close(nc.clim)
nc_close(nc.max)
nc_close(nc.min)

####

pr.clim <- tasmax.clim <- tasmin.clim <- wind.clim <- list()
pr.max <- tasmax.max <- tasmin.max <- wind.max <- list()
pr.min <- tasmax.min <- tasmin.min <- wind.min <- list()
for(i in 1:12){
    pr.clim[[i]] <- pr.clim.mat[,,i]
    tasmax.clim[[i]] <- tasmax.clim.mat[,,i]
    tasmin.clim[[i]] <- tasmin.clim.mat[,,i]
    wind.clim[[i]] <- wind.clim.mat[,,i]
    pr.max[[i]] <- pr.max.mat[,,i]
    tasmax.max[[i]] <- tasmax.max.mat[,,i]
    tasmin.max[[i]] <- tasmin.max.mat[,,i]
    wind.max[[i]] <- wind.max.mat[,,i]
    pr.min[[i]] <- pr.min.mat[,,i]
    tasmax.min[[i]] <- tasmax.min.mat[,,i]
    tasmin.min[[i]] <- tasmin.min.mat[,,i]
    wind.min[[i]] <- wind.min.mat[,,i]
}

####

save(pr.clim, file=paste(output.dir, 'pr.obs.clim', output.suffix,
     '.RData', sep=''))
save(tasmax.clim, file=paste(output.dir, 'tasmax.obs.clim', output.suffix,
     '.RData', sep=''))
save(tasmin.clim, file=paste(output.dir, 'tasmin.obs.clim', output.suffix,
     '.RData', sep=''))
save(wind.clim, file=paste(output.dir, 'wind.obs.clim', output.suffix,
     '.RData', sep=''))

save(pr.max, file=paste(output.dir, 'pr.obs.max', output.suffix,
     '.RData', sep=''))
save(tasmax.max, file=paste(output.dir, 'tasmax.obs.max', output.suffix,
     '.RData', sep=''))
save(tasmin.max, file=paste(output.dir, 'tasmin.obs.max', output.suffix,
     '.RData', sep=''))
save(wind.max, file=paste(output.dir, 'wind.obs.max', output.suffix,
     '.RData', sep=''))

save(pr.min, file=paste(output.dir, 'pr.obs.min', output.suffix,
     '.RData', sep=''))
save(tasmax.min, file=paste(output.dir, 'tasmax.obs.min', output.suffix,
     '.RData', sep=''))
save(tasmin.min, file=paste(output.dir, 'tasmin.obs.min', output.suffix,
     '.RData', sep=''))
save(wind.min, file=paste(output.dir, 'wind.obs.min', output.suffix,
     '.RData', sep=''))

##******************************************************************************
