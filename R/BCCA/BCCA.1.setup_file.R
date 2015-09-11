##******************************************************************************
# Bias Corrected Constructed Analogue (BCCA) downscaling algorithm
# Alex Cannon (acannon@uvic.ca)
# *All wind related lines of code - commented out by Arelia Werner (wernera@uvic.ca)
##******************************************************************************
# Read fine-scale grid and spatially aggregate to GCM grid
##******************************************************************************

ptm <- proc.time()

library(ncdf4)

code.dir <- '/home/ssobie/stat.downscaling/code/QPQM/BCCA/'
source(paste(code.dir,'netcdf.calendar.R',sep=''))

args <- commandArgs(trailingOnly=TRUE)
for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
}

config <- paste(code.dir,'BCCA.set.config',sep='')
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

save(obs.lons, file=paste(output.dir, 'obs.lons', output.suffix,
     '.RData', sep=''))
save(obs.lats, file=paste(output.dir, 'obs.lats', output.suffix,
     '.RData', sep=''))
save(obs.time, file=paste(output.dir, 'obs.time', output.suffix,
     '.RData', sep=''))

save(gcm.lons, file=paste(output.dir, 'gcm.lons', output.suffix,
     '.RData', sep=''))
save(gcm.lats, file=paste(output.dir, 'gcm.lats', output.suffix,
     '.RData', sep=''))

################################################################################
# Figure out which GCM grid boxes are associated with each fine-scale grid point
# Confine search to 10 deg. x 10 deg. neighbourhood

dxy <- 10
mdist <- function(x, y)
    apply(abs(sweep(data.matrix(y), 2, data.matrix(x), '-')), 1, sum)
nn <- list()
for (i in seq_along(obs.lons)) {
    if((i %% 500)==0) cat(i, '')
    gcm.lims <- ((gcm.lons.lats[,1] >= (obs.lons[i]-dxy)) &
                 (gcm.lons.lats[,1] <= (obs.lons[i]+dxy))) &
                ((gcm.lons.lats[,2] >= (obs.lats[i]-dxy)) &
                 (gcm.lons.lats[,2] <= (obs.lats[i]+dxy)))
    gcm.lims <- which(gcm.lims)
    nn.min <- which.min(mdist(c(obs.lons[i], obs.lats[i]),
                        gcm.lons.lats[gcm.lims,]))
    nn[[i]] <- gcm.lims[nn.min]
}
nn <- unlist(nn)

save(nn, file=paste(output.dir, 'nn.index', output.suffix,
                      '.RData', sep=''))

################################################################################

save(pr.gcm.time, file=paste(output.dir, 'pr.gcm.time', output.suffix,
                      '.RData', sep=''))
save(tasmax.gcm.time, file=paste(output.dir, 'tasmax.gcm.time', output.suffix,
                      '.RData', sep=''))
save(tasmin.gcm.time, file=paste(output.dir, 'tasmin.gcm.time', output.suffix,
                      '.RData', sep=''))
save(pr.raw.time, file=paste(output.dir, 'pr.raw.time', output.suffix,
                      '.RData', sep=''))
save(tasmax.raw.time, file=paste(output.dir, 'tasmax.raw.time', output.suffix,
                      '.RData', sep=''))
save(tasmin.raw.time, file=paste(output.dir, 'tasmin.raw.time', output.suffix,
                      '.RData', sep=''))
##******************************************************************************

print('Elapsed Time')
print(proc.time() - ptm)