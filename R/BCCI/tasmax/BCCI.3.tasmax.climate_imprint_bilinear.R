##******************************************************************************
# Bias Corrected Climate Imprint (BCCI) downscaling algorithm
# Alex Cannon (acannon@uvic.ca)
# *All wind related lines of code - commented out by Arelia Werner (wernera@uvic.ca)
##******************************************************************************
# Bilinearly interpolate daily GCM anomalies to fine-scale grid and
# superimpose back onto fine-scale monthly climatologies
##******************************************************************************
##Modified to only downscale tasmax

ptm <- proc.time()

code.dir <- '/home/ssobie/stat.downscaling/code/QPQM/BCCI/'

library(fields)
library(ncdf4)

args <- commandArgs(trailingOnly=TRUE)
for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
}

config <- paste(code.dir,'BCCI.set.config',sep='')
print(readLines(config))
source(config)

  nc <- nc_open(paste(output.dir, 'imp.', output.file, output.suffix, '_tasmax.nc',
                      sep=''), write=TRUE)
  
  tasmax.missing_value <- ncatt_get(nc, varid='tasmax',
                                    attname='missing_value')$value
  tasmax.scale_factor <- ncatt_get(nc, varid='tasmax',
                                   attname='scale_factor')$value
  tasmax.add_offset <- ncatt_get(nc, varid='tasmax',
                                 attname='add_offset')$value

  load(file=paste(output.dir, 'tasmax.dates.gcm', output.suffix, '.RData', sep=''))
  load(file=paste(output.dir, 'lonlat.gcm', output.suffix, '.RData', sep=''))
  load(file=paste(obs.dir, 'lonlat.obs.clim1951-2005.RData', sep=''))
  load(file=paste(output.dir, 'tasmax.raw.time', output.suffix, '.RData', sep=''))

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

  load(file=paste(output.dir, 'tasmax.gcm', output.suffix, '.RData', sep=''))
  load(file=paste(obs.dir, 'tasmax.obs.clim1951-2005.RData', sep=''))
  obs.dim <- dim(tasmax.clim[[1]])

  for(i in seq_along(dates[,1])){
    cat(dates[i,], '\n')
    mn <- dates[i,2]
    ncvar_put(nc=nc, varid='time', vals=raw.time[i], start=i, count=1)

    tasmax.grid <- list(x=lonlat$lon-360, y=lonlat$lat, z=tasmax[,,i])
    bilin.tasmax <- interp.surface(tasmax.grid, obs.coords)
    bilin.tasmax <- bilin.tasmax + c(tasmax.clim[[mn]])
    bilin.tasmax <- round((bilin.tasmax - tasmax.add_offset)/
                          tasmax.scale_factor)
    bilin.tasmax[is.na(bilin.tasmax)] <- tasmax.missing_value
    dim(bilin.tasmax) <- obs.dim
    ncvar_put(nc=nc, varid='tasmax', vals=bilin.tasmax, start=c(1, 1, i),
              count=c(obs.dim[1], obs.dim[2], 1))
    nc_sync(nc)
  }
  rm(tasmax)
  rm(tasmax.clim)
  rm(bilin.tasmax)
  rm(tasmax.grid)
  gc()
##******************************************************************************

print('Elapsed time')
print(proc.time() - ptm)

rm(list=ls())



