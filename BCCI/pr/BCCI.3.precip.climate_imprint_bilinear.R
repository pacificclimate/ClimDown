##******************************************************************************
# Bias Corrected Climate Imprint (BCCI) downscaling algorithm
# Alex Cannon (acannon@uvic.ca)
# *All wind related lines of code - commented out by Arelia Werner (wernera@uvic.ca)
##******************************************************************************
# Bilinearly interpolate daily GCM anomalies to fine-scale grid and
# superimpose back onto fine-scale monthly climatologies
##******************************************************************************
##Modified to analyze only precipitation

ptm <- proc.time()

library(fields)
library(ncdf4)

code.dir <- '/home/ssobie/stat.downscaling/code/QPQM/BCCI/'

config <- paste(code.dir,'BCCI.set.config',sep='')
print(readLines(config))
source(config)

  nc <- nc_open(paste(output.dir, 'imp.', output.file, output.suffix, '_precip.nc',
                      sep=''), write=TRUE)
  
  pr.missing_value <- ncatt_get(nc, varid='pr',
                                attname='missing_value')$value
  pr.scale_factor <- ncatt_get(nc, varid='pr',
                               attname='scale_factor')$value
  pr.add_offset <- ncatt_get(nc, varid='pr',
                             attname='add_offset')$value

  load(file=paste(output.dir, 'pr.dates.gcm', output.suffix, '.RData', sep=''))
  load(file=paste(output.dir, 'lonlat.gcm', output.suffix, '.RData', sep=''))
  load(file=paste(obs.dir, 'lonlat.obs.clim1951-2005.RData', sep=''))
  load(file=paste(output.dir, 'pr.raw.time', output.suffix, '.RData', sep=''))

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

  load(file=paste(output.dir, 'pr.gcm', output.suffix, '.RData', sep=''))
  load(file=paste(obs.dir, 'pr.obs.clim1951-2005.RData', sep=''))
  obs.dim <- dim(pr.clim[[1]])

  for(i in seq_along(dates[,1])){
    cat(dates[i,], '\n')
    mn <- dates[i,2]
    ncvar_put(nc=nc, varid='time', vals=raw.time[i], start=i, count=1)
    ##
    pr.grid <- list(x=lonlat$lon-360, y=lonlat$lat, z=pr[,,i])
    bilin.pr <- interp.surface(pr.grid, obs.coords)
    bilin.pr <- bilin.pr*c(pr.clim[[mn]])
    bilin.pr[bilin.pr > pr.overflow] <- pr.overflow
    bilin.pr <- round((bilin.pr - pr.add_offset)/pr.scale_factor)
    bilin.pr[is.na(bilin.pr)] <- pr.missing_value
    dim(bilin.pr) <- obs.dim
    ncvar_put(nc=nc, varid='pr', vals=bilin.pr, start=c(1, 1, i),
              count=c(obs.dim[1], obs.dim[2], 1))
    nc_sync(nc)
  }
  rm(pr)
  rm(pr.clim)
  rm(bilin.pr)
  rm(pr.grid)
  gc()
##******************************************************************************

print('Elapsed time')
print(proc.time() - ptm)

rm(list=ls())



