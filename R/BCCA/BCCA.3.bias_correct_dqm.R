##******************************************************************************
# Bias Corrected Constructed Analogue (BCCA) downscaling algorithm
# Alex Cannon (acannon@uvic.ca)
# *All wind related lines of code - commented out by Arelia Werner (wernera@uvic.ca)
##******************************************************************************
# Bias correct daily GCM values using detrended quantile mapping algorithm
##******************************************************************************

library(ncdf4)

ptm <- proc.time()
print('Starting step 3')

source('../DQM.R')

args <- commandArgs(trailingOnly=TRUE)
for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
}

config <- paste(code.dir,'BCCA.set.config',sep='')
print(readLines(config))
source(config)

load(paste(output.dir, 'gcm.lons', output.suffix, '.RData', sep=''))
load(paste(output.dir, 'gcm.lats', output.suffix, '.RData', sep=''))

load(paste(output.dir, 'obs.time', output.suffix, '.RData', sep=''))

load(paste(output.dir, 'pr.gcm.time', output.suffix, '.RData', sep=''))
load(paste(output.dir, 'tasmax.gcm.time', output.suffix, '.RData', sep=''))
load(paste(output.dir, 'tasmin.gcm.time', output.suffix, '.RData', sep=''))

h.obs.indices <- obs.time[,1] %in% bc.years
h.gcm.indices <- pr.gcm.time[,1] %in% bc.years
f.gcm.indices <- pr.gcm.time[,1] > max(bc.years)
p.gcm.indices <- pr.gcm.time[,1] < min(bc.years)

##******************************************************************************

# Mask GCM grid points based on spatially aggregated observed land mask
load(paste(output.dir, 'pr.aggregate.one', output.suffix, '.RData', sep=''))
na.mask <- pr.aggregate.one*0 + 1
points <- which(!is.na(na.mask))

##Precipitation
# Read in GCM data
load(paste(output.dir, 'pr.aggregate', output.suffix, '.RData', sep=''))
pr.nc <- nc_open(pr.nc.file)
print(pr.nc)
pr.gcm <- ncvar_get(pr.nc, pr.var)*86400
pr.gcm <- round(pr.gcm, 3)
nc_close(pr.nc)

pr.gcm <- aperm(pr.gcm, c(3, 1, 2))
dim(pr.gcm) <- c(dim(pr.gcm)[1], prod(dim(pr.gcm)[2:3]))
pr.gcm <- sweep(pr.gcm, 2, na.mask, '*')

# Bias correct daily GCM values using detrended quantile mapping

pr.gcm.tmp <- list()
print('PR step')

for(i in seq_along(points)) {
#print(i)
  point <- points[i]
#  cat(point, '\n')
#  print(length(pr.gcm[,point]))
  if (sum(is.na(pr.gcm[,point])) == length(pr.gcm[,point])) {
     pr.gcm.tmp[[i]] <- rep(NA,dim(pr.gcm)[1])
  } else { 
     dqm.pr.tmp <- mnDQM(obs.h=pr.aggregate[h.obs.indices,point],
                      gcm.h=pr.gcm[h.gcm.indices,point],
                      gcm.f=pr.gcm[f.gcm.indices,point],
                      months.obs.h=obs.time[h.obs.indices,2],
                      months.gcm.h=pr.gcm.time[h.gcm.indices,2],
                      months.gcm.f=pr.gcm.time[f.gcm.indices,2],
                      gcm.p=pr.gcm[p.gcm.indices,point],
                      months.gcm.p=pr.gcm.time[p.gcm.indices,2],
                      ratio=TRUE, detrend=detrend.pr, n.max=NULL)
     pr.gcm.tmp[[i]] <- c(dqm.pr.tmp$g.p.bc, dqm.pr.tmp$g.h.bc, dqm.pr.tmp$g.f.bc)
   }
   
}

for(i in seq_along(points)){
  point <- points[i]
  pr.gcm[,point] <- pr.gcm.tmp[[i]]
}

save(pr.gcm, file=paste(output.dir, 'pr.gcm_bc', output.suffix,
     '.RData', sep=''))

rm(pr.gcm)
rm(pr.gcm.tmp)
rm(dqm.pr.tmp)
rm(pr.aggregate)

################################################################################
if (1==1) {
# TASMAX
load(paste(output.dir, 'tasmax.aggregate', output.suffix, '.RData', sep=''))

h.gcm.indices <- tasmax.gcm.time[,1] %in% bc.years
f.gcm.indices <- tasmax.gcm.time[,1] > max(bc.years)
p.gcm.indices <- tasmax.gcm.time[,1] < min(bc.years)

# Read in GCM Data
tasmax.nc <- nc_open(tasmax.nc.file)
tasmax.gcm <- ncvar_get(tasmax.nc, tasmax.var)-273.15
nc_close(tasmax.nc)
##******************************************************************************
tasmax.gcm <- aperm(tasmax.gcm, c(3, 1, 2))
dim(tasmax.gcm) <- c(dim(tasmax.gcm)[1], prod(dim(tasmax.gcm)[2:3]))
tasmax.gcm <- sweep(tasmax.gcm, 2, na.mask, '*')

print('TASMAX step')
##browser()

tasmax.gcm.tmp <- list()

for (i in seq_along(points)) {
  print(i)
  point <- points[i]
  cat(point, '\n')

  if (sum(is.na(tasmax.gcm[,point])) == length(tasmax.gcm[,point])) {
     tasmax.gcm.tmp[[i]] <- rep(NA,dim(tasmax.gcm)[1])
  } else { 
      dqm.tasmax.tmp <- mnDQM(obs.h=tasmax.aggregate[h.obs.indices,point],
                          gcm.h=tasmax.gcm[h.gcm.indices,point],
                          gcm.f=tasmax.gcm[f.gcm.indices,point],
                          months.obs.h=obs.time[h.obs.indices,2],
                          months.gcm.h=tasmax.gcm.time[h.gcm.indices,2],
                          months.gcm.f=tasmax.gcm.time[f.gcm.indices,2],
                          gcm.p=tasmax.gcm[p.gcm.indices,point],
                          months.gcm.p=tasmax.gcm.time[p.gcm.indices,2],
                          ratio=FALSE, detrend=detrend.tx, n.max=NULL)
      tasmax.gcm.tmp[[i]] <- c(dqm.tasmax.tmp$g.p.bc, dqm.tasmax.tmp$g.h.bc, dqm.tasmax.tmp$g.f.bc)
   }
}
for(i in seq_along(points)){
  point <- points[i]
  tasmax.gcm[,point] <- tasmax.gcm.tmp[[i]]
}

save(tasmax.gcm, file=paste(output.dir, 'tasmax.gcm_bc', output.suffix,
     '.RData', sep=''))

rm(tasmax.gcm)
rm(tasmax.gcm.tmp)
rm(dqm.tasmax.tmp)
rm(tasmax.aggregate)

################################################################################
# TASMIN
load(paste(output.dir, 'tasmin.aggregate', output.suffix, '.RData', sep=''))

h.gcm.indices <- tasmin.gcm.time[,1] %in% bc.years
f.gcm.indices <- tasmin.gcm.time[,1] > max(bc.years)
p.gcm.indices <- tasmin.gcm.time[,1] < min(bc.years)

# Read in GCM Data
tasmin.nc <- nc_open(tasmin.nc.file)
tasmin.gcm <- ncvar_get(tasmin.nc, tasmin.var)-273.15
nc_close(tasmin.nc)
##******************************************************************************

tasmin.gcm <- aperm(tasmin.gcm, c(3, 1, 2))
dim(tasmin.gcm) <- c(dim(tasmin.gcm)[1], prod(dim(tasmin.gcm)[2:3]))
tasmin.gcm <- sweep(tasmin.gcm, 2, na.mask, '*')

################################################################################
# Bias correct daily GCM values using detrended quantile mapping

print('TASMIN step')
tasmin.gcm.tmp <- list()

for (i in seq_along(points)) {
  print(i)
  point <- points[i]
  cat(point, '\n')

  if (sum(is.na(tasmin.gcm[,point])) == length(tasmin.gcm[,point])) {
     tasmin.gcm.tmp[[i]] <- rep(NA,dim(tasmin.gcm)[1])
  } else { 
     dqm.tasmin.tmp <- mnDQM(obs.h=tasmin.aggregate[h.obs.indices,point],
                          gcm.h=tasmin.gcm[h.gcm.indices,point],
                          gcm.f=tasmin.gcm[f.gcm.indices,point],
                          months.obs.h=obs.time[h.obs.indices,2],
                          months.gcm.h=tasmin.gcm.time[h.gcm.indices,2],
                          months.gcm.f=tasmin.gcm.time[f.gcm.indices,2],
                          gcm.p=tasmin.gcm[p.gcm.indices,point],
                          months.gcm.p=tasmin.gcm.time[p.gcm.indices,2],
                          ratio=FALSE, detrend=detrend.tn, n.max=NULL)
     tasmin.gcm.tmp[[i]] <- c(dqm.tasmin.tmp$g.p.bc, dqm.tasmin.tmp$g.h.bc, dqm.tasmin.tmp$g.f.bc)
  }
}
print('END step')

for(i in seq_along(points)){
  point <- points[i]
  tasmin.gcm[,point] <- tasmin.gcm.tmp[[i]]
}

save(tasmin.gcm, file=paste(output.dir, 'tasmin.gcm_bc', output.suffix,
     '.RData', sep=''))

rm(tasmin.gcm)
rm(tasmin.gcm.tmp)
rm(dqm.tasmin.tmp)
rm(tasmin.aggregate)



################################################################################


}

print('Elapsed time')
print(proc.time()-ptm)
