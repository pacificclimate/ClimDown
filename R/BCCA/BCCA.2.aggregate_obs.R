##******************************************************************************
# Bias Corrected Constructed Analogues (BCCA) downscaling algorithm
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

##******************************************************************************
# Read fine-scale and GCM grid dimensions

load(file=paste(output.dir, 'gcm.lons', output.suffix,
                      '.RData', sep=''))
load(file=paste(output.dir, 'obs.time', output.suffix,
                      '.RData', sep=''))
load(file=paste(output.dir, 'nn.index', output.suffix,
                      '.RData', sep=''))

gridpoints <- sort(unique(nn))
cat('\n')

obs.lon <- ncvar_get(nc.obs, 'lon')
obs.lat <- ncvar_get(nc.obs, 'lat')
n.lon <- length(obs.lon)
n.lat <- length(obs.lat)

################################################################################
# Spatially aggregate the fine-scale data to the GCM grid

pr.aggregate <- tasmax.aggregate <- tasmin.aggregate <-
      matrix(NA, nrow=nrow(obs.time), ncol=length(gcm.lons))
        #  wind.aggregate <-

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


for(i in seq_along(i.starts)){
    cat(obs.time[i.starts[i],], '\n')
    pr.obs <- ncvar_get(nc.obs, varid='pr', start=c(1, 1, i.starts[i]),
                                                   count=c(n.lon, n.lat, i.lengths[i]))
    dim(pr.obs) <- c(prod(dim(pr.obs)[1:2]), dim(pr.obs)[3])
    pr.agg <- matrix(NA, nrow=i.lengths[i], ncol=length(gcm.lons))
    pr.all.agg <- all.agg.fxn(gridpoints,nn,pr.obs)
    pr.agg[,gridpoints] <- pr.all.agg
    pr.aggregate[i.starts[i]:(i.starts[i]+i.lengths[i]-1),] <- pr.agg
}

save(pr.aggregate, file=paste(output.dir, 'pr.aggregate', output.suffix,
                          '.RData', sep=''))
pr.aggregate.one <- pr.aggregate[1,]
save(pr.aggregate.one, file=paste(output.dir, 'pr.aggregate.one', output.suffix,
                          '.RData', sep=''))
rm(pr.obs)
rm(pr.agg)
rm(pr.all.agg)
rm(pr.aggregate)
gc()

for(i in seq_along(i.starts)){
    cat(obs.time[i.starts[i],], '\n')
    tasmax.obs <- ncvar_get(nc.obs, varid='tasmax', start=c(1, 1, i.starts[i]),
                                                       count=c(n.lon, n.lat, i.lengths[i]))
    dim(tasmax.obs) <- c(prod(dim(tasmax.obs)[1:2]), dim(tasmax.obs)[3])
    tasmax.agg <- matrix(NA, nrow=i.lengths[i], ncol=length(gcm.lons))
    tasmax.all.agg <- all.agg.fxn(gridpoints,nn,tasmax.obs)
    tasmax.agg[,gridpoints] <- tasmax.all.agg
    tasmax.aggregate[i.starts[i]:(i.starts[i]+i.lengths[i]-1),] <- tasmax.agg
}
save(tasmax.aggregate, file=paste(output.dir, 'tasmax.aggregate', output.suffix,
                              '.RData', sep=''))
tasmax.aggregate.one <- tasmax.aggregate[1,]
save(tasmax.aggregate.one, file=paste(output.dir, 'tasmax.aggregate.one', output.suffix,
                              '.RData', sep=''))

rm(tasmax.obs)
rm(tasmax.agg)
rm(tasmax.all.agg)
rm(tasmax.aggregate)
gc()

for(i in seq_along(i.starts)){
    cat(obs.time[i.starts[i],], '\n')
    tasmin.obs <- ncvar_get(nc.obs, varid='tasmin', start=c(1, 1, i.starts[i]),
                                                       count=c(n.lon, n.lat, i.lengths[i]))
    dim(tasmin.obs) <- c(prod(dim(tasmin.obs)[1:2]), dim(tasmin.obs)[3])
    tasmin.agg <- matrix(NA, nrow=i.lengths[i], ncol=length(gcm.lons))
    tasmin.all.agg <- all.agg.fxn(gridpoints,nn,tasmin.obs)
    tasmin.agg[,gridpoints] <- tasmin.all.agg
    tasmin.aggregate[i.starts[i]:(i.starts[i]+i.lengths[i]-1),] <- tasmin.agg
}
save(tasmin.aggregate, file=paste(output.dir, 'tasmin.aggregate', output.suffix,
                              '.RData', sep=''))
tasmin.aggregate.one <- tasmin.aggregate[1,]
save(tasmin.aggregate.one, file=paste(output.dir, 'tasmin.aggregate.one', output.suffix,
                              '.RData', sep=''))
rm(tasmin.obs)
rm(tasmin.agg)
rm(tasmin.all.agg)
rm(tasmin.aggregate)
gc()

print('Elapsed time')
print(proc.time() - ptm)