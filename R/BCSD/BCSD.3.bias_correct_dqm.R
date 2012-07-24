##******************************************************************************
# Bias Corrected Spatial Disaggregation (BCSD) downscaling algorithm
# Alex Cannon (acannon@uvic.ca)
##******************************************************************************
# Bias correct monthly GCM values using detrended quantile mapping algorithm
##******************************************************************************

rm(list=ls())
source('DQM.R')
library(ncdf4)
library(RNetCDF)
library(doMC)

config <- commandArgs(trailingOnly=TRUE)
if(length(config)==0) config <- 'BCSD.config'
print(readLines(config))
source(config)

registerDoMC(cores=mc.cores)

load(paste(output.dir, 'gcm.lons', output.suffix, '.RData', sep=''))
load(paste(output.dir, 'gcm.lats', output.suffix, '.RData', sep=''))
load(paste(output.dir, 'gcm.time', output.suffix, '.RData', sep=''))
load(paste(output.dir, 'obs.time', output.suffix, '.RData', sep=''))
load(paste(output.dir, 'pr.aggregate', output.suffix, '.RData', sep=''))
load(paste(output.dir, 'tasmin.aggregate', output.suffix, '.RData', sep=''))
load(paste(output.dir, 'tasmax.aggregate', output.suffix, '.RData', sep=''))
load(paste(output.dir, 'wind.aggregate', output.suffix, '.RData', sep=''))

pr.nc <- nc_open(pr.nc.file)
tasmax.nc <- nc_open(tasmax.nc.file)
tasmin.nc <- nc_open(tasmin.nc.file)
uwnd.nc <- nc_open(uwnd.nc.file)
vwnd.nc <- nc_open(vwnd.nc.file)

##******************************************************************************
# Read in GCM data

pr.gcm <- ncvar_get(pr.nc, pr.var)*86400
pr.gcm <- round(pr.gcm, 3)
tasmax.gcm <- ncvar_get(tasmax.nc, tasmax.var)-273.15
tasmin.gcm <- ncvar_get(tasmin.nc, tasmin.var)-273.15
uwnd.gcm <- ncvar_get(uwnd.nc, uwnd.var)
vwnd.gcm <- ncvar_get(vwnd.nc, vwnd.var)
wind.gcm <- sqrt(uwnd.gcm^2 + vwnd.gcm^2)
rm(uwnd.gcm, vwnd.gcm)

nc_close(pr.nc)
nc_close(tasmax.nc)
nc_close(tasmin.nc)
nc_close(uwnd.nc)
nc_close(vwnd.nc)

################################################################################
# Mask GCM grid points based on spatially aggregated observed land mask

na.mask <- pr.aggregate[1,]*0 + 1

pr.gcm <- aperm(pr.gcm, c(3, 1, 2))
dim(pr.gcm) <- c(dim(pr.gcm)[1], prod(dim(pr.gcm)[2:3]))
pr.gcm <- sweep(pr.gcm, 2, na.mask, '*')

tasmax.gcm <- aperm(tasmax.gcm, c(3, 1, 2))
dim(tasmax.gcm) <- c(dim(tasmax.gcm)[1], prod(dim(tasmax.gcm)[2:3]))
tasmax.gcm <- sweep(tasmax.gcm, 2, na.mask, '*')

tasmin.gcm <- aperm(tasmin.gcm, c(3, 1, 2))
dim(tasmin.gcm) <- c(dim(tasmin.gcm)[1], prod(dim(tasmin.gcm)[2:3]))
tasmin.gcm <- sweep(tasmin.gcm, 2, na.mask, '*')

wind.gcm <- aperm(wind.gcm, c(3, 1, 2))
dim(wind.gcm) <- c(dim(wind.gcm)[1], prod(dim(wind.gcm)[2:3]))
wind.gcm <- sweep(wind.gcm, 2, na.mask, '*')

################################################################################
# Calculate monthly means for aggregated fine-scale data

cat('Calculating monthly means for aggregated fine-scale data...\n')

yr.obs <- as.character(obs.time[,1])
mn.obs <- as.character(obs.time[,2])
for(i in seq_along(mn.obs)) if(nchar(mn.obs[i])==1)
    mn.obs[i] <- paste('0', mn.obs[i], sep='')
yr.mn.obs <- paste(yr.obs, mn.obs, sep='-')

pr.aggregate.mn <- do.call(rbind, mclapply(split(as.data.frame(
                           pr.aggregate), yr.mn.obs), colMeans))
tasmax.aggregate.mn <- do.call(rbind, mclapply(split(as.data.frame(
                               tasmax.aggregate), yr.mn.obs), colMeans))
tasmin.aggregate.mn <- do.call(rbind, mclapply(split(as.data.frame(
                               tasmin.aggregate), yr.mn.obs), colMeans))
wind.aggregate.mn <- do.call(rbind, mclapply(split(as.data.frame(
                             wind.aggregate), yr.mn.obs), colMeans))

################################################################################
# Calculate monthly means for GCM data

yr.gcm <- as.character(gcm.time[,1])
mn.gcm <- as.character(gcm.time[,2])
for(i in seq_along(mn.gcm)) if(nchar(mn.gcm[i])==1)
    mn.gcm[i] <- paste('0', mn.gcm[i], sep='')
yr.mn.gcm <- paste(yr.gcm, mn.gcm, sep='-')

pr.gcm.mn <- do.call(rbind, mclapply(split(as.data.frame(pr.gcm),
                     yr.mn.gcm), colMeans))
tasmax.gcm.mn <- do.call(rbind, mclapply(split(as.data.frame(tasmax.gcm),
                         yr.mn.gcm), colMeans))
tasmin.gcm.mn <- do.call(rbind, mclapply(split(as.data.frame(tasmin.gcm),
                         yr.mn.gcm), colMeans))
wind.gcm.mn <- do.call(rbind, mclapply(split(as.data.frame(wind.gcm),
                       yr.mn.gcm), colMeans))

gcm.time.mn <- do.call(rbind, mclapply(split(as.data.frame(gcm.time),
                       yr.mn.gcm), colMeans))
obs.time.mn <- do.call(rbind, mclapply(split(as.data.frame(obs.time),
                       yr.mn.obs), colMeans))

h.obs.indices <- obs.time.mn[,1] %in% bc.years
h.gcm.indices <- gcm.time.mn[,1] %in% bc.years
f.gcm.indices <- gcm.time.mn[,1] > max(bc.years)
p.gcm.indices <- gcm.time.mn[,1] < min(bc.years)

################################################################################
# Bias correct monthly GCM values using detrended quantile mapping

points <- which(!is.na(na.mask))
dqm.all <- foreach(i = seq_along(points)) %dopar% {
    point <- points[i]
    cat(point, '\n')
    dqm.pr <- mnDQM(obs.h=pr.aggregate.mn[h.obs.indices,point],
                    gcm.h=pr.gcm.mn[h.gcm.indices,point],
                    gcm.f=pr.gcm.mn[f.gcm.indices,point],
                    months.obs.h=obs.time.mn[h.obs.indices,2],
                    months.gcm.h=gcm.time.mn[h.gcm.indices,2],
                    months.gcm.f=gcm.time.mn[f.gcm.indices,2],
                    gcm.p=pr.gcm.mn[p.gcm.indices,point],
                    months.gcm.p=gcm.time.mn[p.gcm.indices,2],
                    ratio=TRUE, detrend=detrend.pr, n.max=NULL)
    pr.gcm.mn.tmp <- c(dqm.pr$g.p.bc, dqm.pr$g.h.bc, dqm.pr$g.f.bc)
    dqm.tasmax <- mnDQM(obs.h=tasmax.aggregate.mn[h.obs.indices,point],
                        gcm.h=tasmax.gcm.mn[h.gcm.indices,point],
                        gcm.f=tasmax.gcm.mn[f.gcm.indices,point],
                        months.obs.h=obs.time.mn[h.obs.indices,2],
                        months.gcm.h=gcm.time.mn[h.gcm.indices,2],
                        months.gcm.f=gcm.time.mn[f.gcm.indices,2],
                        gcm.p=tasmax.gcm.mn[p.gcm.indices,point],
                        months.gcm.p=gcm.time.mn[p.gcm.indices,2],
                        ratio=FALSE, detrend=detrend.tx, n.max=NULL)
    tasmax.gcm.mn.tmp <- c(dqm.tasmax$g.p.bc, dqm.tasmax$g.h.bc,
                           dqm.tasmax$g.f.bc)
    dqm.tasmin <- mnDQM(obs.h=tasmin.aggregate.mn[h.obs.indices,point],
                        gcm.h=tasmin.gcm.mn[h.gcm.indices,point],
                        gcm.f=tasmin.gcm.mn[f.gcm.indices,point],
                        months.obs.h=obs.time.mn[h.obs.indices,2],
                        months.gcm.h=gcm.time.mn[h.gcm.indices,2],
                        months.gcm.f=gcm.time.mn[f.gcm.indices,2],
                        gcm.p=tasmin.gcm.mn[p.gcm.indices,point],
                        months.gcm.p=gcm.time.mn[p.gcm.indices,2],
                        ratio=FALSE, detrend=detrend.tn, n.max=NULL)
    tasmin.gcm.mn.tmp <- c(dqm.tasmin$g.p.bc, dqm.tasmin$g.h.bc,
                           dqm.tasmin$g.f.bc)
    dqm.wind <- mnDQM(obs.h=wind.aggregate.mn[h.obs.indices,point],
                      gcm.h=wind.gcm.mn[h.gcm.indices,point],
                      gcm.f=wind.gcm.mn[f.gcm.indices,point],
                      months.obs.h=obs.time.mn[h.obs.indices,2],
                      months.gcm.h=gcm.time.mn[h.gcm.indices,2],
                      months.gcm.f=gcm.time.mn[f.gcm.indices,2],
                      gcm.p=wind.gcm.mn[p.gcm.indices,point],
                      months.gcm.p=gcm.time.mn[p.gcm.indices,2],
                      ratio=TRUE, detrend=detrend.wind, n.max=NULL)
    wind.gcm.mn.tmp <- c(dqm.wind$g.p.bc, dqm.wind$g.h.bc, dqm.wind$g.f.bc)
    cbind(pr.gcm.mn.tmp, tasmax.gcm.mn.tmp, tasmin.gcm.mn.tmp, wind.gcm.mn.tmp)
}

for(i in seq_along(points)){
    point <- points[i]
    pr.gcm.mn[,point] <- dqm.all[[i]][,1]
    tasmax.gcm.mn[,point] <- dqm.all[[i]][,2]
    tasmin.gcm.mn[,point] <- dqm.all[[i]][,3]
    wind.gcm.mn[,point] <- dqm.all[[i]][,4]
}

################################################################################

save(obs.time.mn, file=paste(output.dir, 'obs.time.mn', output.suffix,
     '.RData', sep=''))
save(pr.aggregate.mn, file=paste(output.dir, 'pr.aggregate.mn', output.suffix,
     '.RData', sep=''))
save(tasmax.aggregate.mn, file=paste(output.dir, 'tasmax.aggregate.mn',
     output.suffix, '.RData', sep=''))
save(tasmin.aggregate.mn, file=paste(output.dir, 'tasmin.aggregate.mn',
     output.suffix, '.RData', sep=''))
save(wind.aggregate.mn, file=paste(output.dir, 'wind.aggregate.mn',
     output.suffix, '.RData', sep=''))

save(gcm.time.mn, file=paste(output.dir, 'gcm.time.mn', output.suffix,
     '.RData', sep=''))
save(pr.gcm.mn, file=paste(output.dir, 'pr.gcm_bc.mn', output.suffix,
     '.RData', sep=''))
save(tasmax.gcm.mn, file=paste(output.dir, 'tasmax.gcm_bc.mn', output.suffix,
     '.RData', sep=''))
save(tasmin.gcm.mn, file=paste(output.dir, 'tasmin.gcm_bc.mn', output.suffix,
     '.RData', sep=''))
save(wind.gcm.mn, file=paste(output.dir, 'wind.gcm_bc.mn', output.suffix,
     '.RData', sep=''))

##******************************************************************************
