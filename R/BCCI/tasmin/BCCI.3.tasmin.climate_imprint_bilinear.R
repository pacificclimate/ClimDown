##******************************************************************************
# Bias Corrected Climate Imprint (BCCI) downscaling algorithm
# Alex Cannon (acannon@uvic.ca)
# *All wind related lines of code - commented out by Arelia Werner (wernera@uvic.ca)
##******************************************************************************
# Bilinearly interpolate daily GCM anomalies to fine-scale grid and
# superimpose back onto fine-scale monthly climatologies
##******************************************************************************

ptm <- proc.time()

library(fields)
library(ncdf4)

code.dir <- '/home/ssobie/stat.downscaling/code/QPQM/BCCI/'

args <- commandArgs(trailingOnly=TRUE)
for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
}

config <- paste(code.dir,'BCCI.set.config',sep='')
print(readLines(config))
source(config)
  
  nc <- nc_open(paste(output.dir, 'imp.', output.file, output.suffix, '_tasmin.nc',
                      sep=''), write=TRUE)

  tasmin.missing_value <- ncatt_get(nc, varid='tasmin',
                                    attname='missing_value')$value
  tasmin.scale_factor <- ncatt_get(nc, varid='tasmin',
                                   attname='scale_factor')$value
  tasmin.add_offset <- ncatt_get(nc, varid='tasmin',
                                 attname='add_offset')$value

  load(file=paste(output.dir, 'tasmin.dates.gcm', output.suffix, '.RData', sep=''))
  load(file=paste(output.dir, 'lonlat.gcm', output.suffix, '.RData', sep=''))
  load(file=paste(obs.dir, 'lonlat.obs.clim1951-2005.RData', sep=''))
  load(file=paste(output.dir, 'tasmin.raw.time', output.suffix, '.RData', sep=''))

##******************************************************************************
# Read observed coordinates
  
  obs.lats <- c(matrix(lonlat.clim$lat, ncol=length(lonlat.clim$lat),
                       nrow=length(lonlat.clim$lon), byrow=TRUE))
  obs.lons <- c(matrix(lonlat.clim$lon, ncol=length(lonlat.clim$lat),
                       nrow=length(lonlat.clim$lon)))
  obs.coords <- cbind(obs.lons, obs.lats)

################################################################################
# Interpolate GCM anomalies to fine-scale grid and superimpose onto fine-scale
# climatologies

  load(file=paste(output.dir, 'tasmin.gcm', output.suffix, '.RData', sep=''))
  load(file=paste(obs.dir, 'tasmin.obs.clim1951-2005.RData', sep=''))
  obs.dim <- dim(tasmin.clim[[1]])

  for(i in seq_along(dates[,1])){
    cat(dates[i,], '\n')
    mn <- dates[i,2]    
    ncvar_put(nc=nc, varid='time', vals=raw.time[i], start=i, count=1)
    tasmin.grid <- list(x=lonlat$lon-360, y=lonlat$lat, z=tasmin[,,i])
    bilin.tasmin <- interp.surface(tasmin.grid, obs.coords)
    bilin.tasmin <- bilin.tasmin + c(tasmin.clim[[mn]])
    bilin.tasmin <- round((bilin.tasmin - tasmin.add_offset)/
                          tasmin.scale_factor)
    bilin.tasmin[is.na(bilin.tasmin)] <- tasmin.missing_value
    dim(bilin.tasmin) <- obs.dim
    ncvar_put(nc=nc, varid='tasmin', vals=bilin.tasmin, start=c(1, 1, i),
              count=c(obs.dim[1], obs.dim[2], 1))

    nc_sync(nc)
  }
  rm(tasmin)
  rm(tasmin.clim)
  rm(bilin.tasmin)
  rm(tasmin.grid)
  gc()
  nc_close(nc)

##******************************************************************************


print('Elapsed time')
print(proc.time - ptm)

rm(list=ls())



