##******************************************************************************
# Bias Corrected Spatial Disaggregation (BCSD) downscaling algorithm
# Alex Cannon (acannon@uvic.ca)
##******************************************************************************
# Read fine-scale grid and spatially aggregate to GCM grid
##******************************************************************************

rm(list=ls())
library(ncdf4)
library(RNetCDF)
library(doMC)
source('netcdf.calendar.R')

config <- commandArgs(trailingOnly=TRUE)
if(length(config)==0) config <- 'BCSD.config'
print(readLines(config))
source(config)

registerDoMC(cores=mc.cores)
nc.obs <- nc_open(nc.obs.file)
nc.gcm <- nc_open(pr.nc.file)

##******************************************************************************
# Read fine-scale and GCM grid dimensions

obs.lon <- ncvar_get(nc.obs, 'lon')
obs.lat <- ncvar_get(nc.obs, 'lat')
n.lon <- length(obs.lon)
n.lat <- length(obs.lat)

obs.lats <- matrix(obs.lat, nrow=n.lon, ncol=n.lat, byrow=TRUE)
obs.lons <- matrix(obs.lon, nrow=n.lon, ncol=n.lat)
obs.time <- netcdf.calendar(nc.obs)

gcm.lon <- ncvar_get(nc.gcm, 'lon')-360
gcm.lat <- ncvar_get(nc.gcm, 'lat')
gcm.lats <- matrix(gcm.lat, ncol=length(gcm.lat), nrow=length(gcm.lon),
                   byrow=TRUE)
gcm.lons <- matrix(gcm.lon, ncol=length(gcm.lat), nrow=length(gcm.lon))
gcm.lons.lats <- cbind(c(gcm.lons), c(gcm.lats))
gcm.time <- netcdf.calendar(nc.gcm)

nc_close(nc.gcm)

################################################################################
# Figure out which GCM grid boxes are associated with each fine-scale grid point
# Confine search to 15 deg. x 15 deg. neighbourhood

dxy <- 15
mdist <- function(x, y)
    apply(abs(sweep(data.matrix(y), 2, data.matrix(x), '-')), 1, sum)
nn <- foreach(i = seq_along(obs.lons)) %dopar% {
    if((i %% 500)==0) cat(i, '')
    gcm.lims <- ((gcm.lons.lats[,1] >= (obs.lons[i]-dxy)) &
                 (gcm.lons.lats[,1] <= (obs.lons[i]+dxy))) &
                ((gcm.lons.lats[,2] >= (obs.lats[i]-dxy)) &
                 (gcm.lons.lats[,2] <= (obs.lats[i]+dxy)))
    gcm.lims <- which(gcm.lims)
    nn.min <- which.min(mdist(c(obs.lons[i], obs.lats[i]),
                        gcm.lons.lats[gcm.lims,]))
    gcm.lims[nn.min]
}
nn <- unlist(nn)
gridpoints <- sort(unique(nn))
cat('\n')

################################################################################
# Spatially aggregate the fine-scale data to the GCM grid

pr.aggregate <- tasmax.aggregate <- tasmin.aggregate <- wind.aggregate <-
    matrix(NA, nrow=nrow(obs.time), ncol=length(gcm.lons))

i.starts <- sapply(split(seq_along(obs.time[,1]), obs.time[,1]), min)
i.lengths <- sapply(split(seq_along(obs.time[,1]), obs.time[,1]), length)

for(i in seq_along(i.starts)){
    cat(obs.time[i.starts[i],], '\n')
    pr.obs <- ncvar_get(nc.obs, varid='pr', start=c(1, 1, i.starts[i]),
                        count=c(n.lon, n.lat, i.lengths[i]))
    tasmax.obs <- ncvar_get(nc.obs, varid='tasmax', start=c(1, 1, i.starts[i]),
                            count=c(n.lon, n.lat, i.lengths[i]))
    tasmin.obs <- ncvar_get(nc.obs, varid='tasmin', start=c(1, 1, i.starts[i]),
                            count=c(n.lon, n.lat, i.lengths[i]))
    wind.obs <- ncvar_get(nc.obs, varid='wind', start=c(1, 1, i.starts[i]),
                          count=c(n.lon, n.lat, i.lengths[i]))
    dim(pr.obs) <- c(prod(dim(pr.obs)[1:2]), dim(pr.obs)[3])
    dim(tasmax.obs) <- c(prod(dim(tasmax.obs)[1:2]), dim(tasmax.obs)[3])
    dim(tasmin.obs) <- c(prod(dim(tasmin.obs)[1:2]), dim(tasmin.obs)[3])
    dim(wind.obs) <- c(prod(dim(wind.obs)[1:2]), dim(wind.obs)[3])
    pr.agg <- tasmax.agg <- tasmin.agg <- wind.agg <- matrix(NA, nrow=i.lengths[i],
                                                             ncol=length(gcm.lons))
    all.agg <- foreach(j=1:length(gridpoints)) %dopar% {
        point <- gridpoints[j]
        cbind(apply(pr.obs[nn==point,], 2, mean, trim=0.1, na.rm=TRUE),
              apply(tasmax.obs[nn==point,], 2, mean, trim=0.1, na.rm=TRUE),
              apply(tasmin.obs[nn==point,], 2, mean, trim=0.1, na.rm=TRUE),
              apply(wind.obs[nn==point,], 2, mean, trim=0.1, na.rm=TRUE))
    }
    all.agg <- do.call(cbind, all.agg)
    pr.all.agg <- all.agg[,c(TRUE, FALSE, FALSE, FALSE)]
    tasmax.all.agg <- all.agg[,c(FALSE, TRUE, FALSE, FALSE)]
    tasmin.all.agg <- all.agg[,c(FALSE, FALSE, TRUE, FALSE)]
    wind.all.agg <- all.agg[,c(FALSE, FALSE, FALSE, TRUE)]
    pr.agg[,gridpoints] <- pr.all.agg
    tasmax.agg[,gridpoints] <- tasmax.all.agg
    tasmin.agg[,gridpoints] <- tasmin.all.agg
    wind.agg[,gridpoints] <- wind.all.agg
    pr.agg[is.nan(pr.agg)] <- NA
    tasmax.agg[is.nan(tasmax.agg)] <- NA
    tasmin.agg[is.nan(tasmin.agg)] <- NA
    wind.agg[is.nan(wind.agg)] <- NA
    pr.aggregate[i.starts[i]:(i.starts[i]+i.lengths[i]-1),] <- pr.agg
    tasmax.aggregate[i.starts[i]:(i.starts[i]+i.lengths[i]-1),] <- tasmax.agg
    tasmin.aggregate[i.starts[i]:(i.starts[i]+i.lengths[i]-1),] <- tasmin.agg
    wind.aggregate[i.starts[i]:(i.starts[i]+i.lengths[i]-1),] <- wind.agg
}
nc_close(nc.obs)

################################################################################

save(gcm.lons, file=paste(output.dir, 'gcm.lons', output.suffix,
     '.RData', sep=''))
save(gcm.lats, file=paste(output.dir, 'gcm.lats', output.suffix,
     '.RData', sep=''))
save(gcm.time, file=paste(output.dir, 'gcm.time', output.suffix,
     '.RData', sep=''))
save(obs.lons, file=paste(output.dir, 'obs.lons', output.suffix,
     '.RData', sep=''))
save(obs.lats, file=paste(output.dir, 'obs.lats', output.suffix,
     '.RData', sep=''))
save(obs.time, file=paste(output.dir, 'obs.time', output.suffix,
     '.RData', sep=''))
save(pr.aggregate, file=paste(output.dir, 'pr.aggregate', output.suffix,
     '.RData', sep=''))
save(tasmax.aggregate, file=paste(output.dir, 'tasmax.aggregate',
     output.suffix, '.RData', sep=''))
save(tasmin.aggregate, file=paste(output.dir, 'tasmin.aggregate',
     output.suffix, '.RData', sep=''))
save(wind.aggregate, file=paste(output.dir, 'wind.aggregate', output.suffix,
     '.RData', sep=''))

##******************************************************************************
