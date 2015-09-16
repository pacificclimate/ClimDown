##******************************************************************************
# Bias Corrected Constructed Analogue (BCCA) downscaling algorithm
# Alex Cannon (acannon@uvic.ca)
# *All wind related lines of code - commented out by Arelia Werner (wernera@uvic.ca)
##******************************************************************************
# Read fine-scale grid and spatially aggregate to GCM grid
##******************************************************************************

library(ncdf4)
source('../netcdf.calendar.R')
source('../bisect.R')

ptm <- proc.time()

print("Starting step 1")

# "--args gcm='${gcm}' rcp='${rcp}' run='${run}'"
args <- commandArgs(trailingOnly=TRUE)
for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
}

config <- 'BCCA.set.config'
print(readLines(config))
source(config)

nc.obs <- nc_open(nc.obs.file)
nc.gcm <- nc_open(pr.nc.file)
nc.tasmax.gcm <- nc_open(tasmax.nc.file)
nc.tasmin.gcm <- nc_open(tasmin.nc.file)

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

pr.gcm.time <- netcdf.calendar(nc.gcm)
tasmax.gcm.time <- netcdf.calendar(nc.tasmax.gcm)
tasmin.gcm.time <- netcdf.calendar(nc.tasmin.gcm)

pr.raw.time <- ncvar_get(nc.gcm,'time')
tasmax.raw.time <- ncvar_get(nc.tasmax.gcm,'time')
tasmin.raw.time <- ncvar_get(nc.tasmin.gcm,'time')

nc_close(nc.gcm)
nc_close(nc.tasmax.gcm)
nc_close(nc.tasmin.gcm)

################################################################################
# Figure out which GCM grid boxes are associated with each fine-scale grid point

grid.mapping <- regrid.coarse.to.fine(gcm.lats, gcm.lons, obs.lats, obs.lons)

print('Step 1 Elapsed Time')
print(proc.time() - ptm)

##******************************************************************************
# Bias Corrected Constructed Analogues (BCCA) downscaling algorithm
# Alex Cannon (acannon@uvic.ca)
# *All wind related lines of code - commented out by Arelia Werner (wernera@uvic.ca)
##******************************************************************************
# Read fine-scale grid and spatially aggregate to GCM grid
##******************************************************************************

print('Starting step 2')
ptm <- proc.time()

#gridpoints <- sort(unique(nn))
cat('\n')

################################################################################
# Spatially aggregate the fine-scale data to the GCM grid

aggregate <- matrix(NA, nrow=nrow(obs.time), ncol=length(gcm.lons))

i.starts <- sapply(split(seq_along(obs.time[,1]), obs.time[,1]), min)
i.lengths <- sapply(split(seq_along(obs.time[,1]), obs.time[,1]), length)


all.agg.fxn <- function(gridpoints,nn,var.obs) {
  all.agg <- matrix(NA,nrow=dim(var.obs)[2],ncol=length(gridpoints))
  for (j in 1:length(gridpoints)) {
    point <- gridpoints[j]
    all.agg[,j] <- apply(var.obs[nn==point,,drop=FALSE], 2, mean, trim=0.1, na.rm=TRUE)
  }
  return(all.agg)
}



for (varid in c('pr', 'tasmax', 'tasmin')) {
    for(i in seq_along(i.starts)){
        cat(obs.time[i.starts[i],], '\n')
        obs <- ncvar_get(nc.obs, varid=varid, start=c(1, 1, i.starts[i]),
                         count=c(n.lon, n.lat, i.lengths[i]))
        dim(obs) <- c(prod(dim(obs)[1:2]), dim(obs)[3])
        agg <- matrix(NA, nrow=i.lengths[i], ncol=length(gcm.lons))
        all.agg <- all.agg.fxn(gridpoints,nn,obs)
        agg[,gridpoints] <- all.agg
        aggregate[i.starts[i]:(i.starts[i]+i.lengths[i]-1),] <- agg
    }

    save(aggregate, file=paste(output.dir, paste(varid, 'aggregate', sep='.'), output.suffix,
                        '.RData', sep=''))
    aggregate.one <- aggregate[1,]
    save(aggregate.one, file=paste(output.dir, paste(varid, 'aggregate.one', sep='.'), output.suffix,
                            '.RData', sep=''))
    rm(obs)
    rm(agg)
    rm(all.agg)
    rm(aggregate)
    gc()
}

print('Elapsed time')
print(proc.time() - ptm)
