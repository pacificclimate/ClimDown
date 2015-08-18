#!/usr/bin/Rscript

################################################################################

rm(list=ls())

ptm <- proc.time()

args <- c('pr',
          '/home/data/climate/downscale/CMIP5/BCCA/pr+tasmax+tasmin_day_BCCA+ANUSPLIN300+CanESM2_historical+rcp45_r1i1p1_19500101-21001231.h5',
          '/home/data/climate/downscale/CMIP5/nobackup/QPQM/done/pr+tasmax+tasmin_day_QPQM+ANUSPLIN300+CanESM2_historical+rcp45_r1i1p1_19500101-21001231.nc',
          '/datasets/projects-rci-scratch/test/stage_4_canesm.nc','2')

## Help section
if('--help' %in% args || length(args)==0) {
    cat('
      Quantile perturbation quantile mapping

      Usage:
      ./pr+tasmax+tasmin_qpqm_netcdf.R [var.name] [obs_file] [gcm_file] [output_file] [cores]

      Arguments:
      var.name - One of pr, tasmax, tasmin
      obs_file - Gridded historical observations
      gcm_file - GCM simulations interpolated to the obs_file grid
      output_file - The file to create (or overwrite) with the bias corrected outputs
      cores - The number of processor cores

      All files should have the same spatial domain.\n\n')

      q(save='no')
}

library(ncdf4)
library(abind)
library(PCICt)
library(Rmpi)
library(doMPI)	

reorder <- function(x,ix) {
  rx <- x[ix]
  return(rx)
}


cl <- startMPIcluster(count=args[5])
registerDoMPI(cl)
exportDoMPI(cl,c('reorder'))
#exportDoMPI(cl, c('QPQM', 'tQPQM'))
mpi.options <- list(info=FALSE, profile=FALSE)
cat('cores =', getDoParWorkers(), '\n')

pr.offset <- 0.
pr.scale <- 1.
tasmax.offset <- 0.
tasmax.scale <- 1.
tasmin.offset <- 0.
tasmin.scale <- 1.

# Read in the input and output files
##var.name <- 'pr'

var.name <- args[1]
bcca.file <- args[2]
qpqm.file <- args[3]
out.file <- args[4]

print('Copying File')
#file.copy(from=qpqm.file,out.file,overwrite=T)

qpqm <- nc_open(qpqm.file)
bcca <- nc_open(bcca.file)
out  <- nc_open(out.file,write=TRUE)

lat <- qpqm$dim$lat$vals
nlat <- length(lat)
lon <- qpqm$dim$lon$vals
nlon <- length(lon)

dates.qpqm <- as.PCICt(strsplit(qpqm$dim$time$units, ' ')[[1]][3],
                      qpqm$dim$time$calendar) + qpqm$dim$time$vals*86400
dates.qpqm <-  apply(do.call(rbind, strsplit(format(dates.qpqm, '%Y %m %d'),
                    ' ')), 2, as.integer)

dates.bcca <- as.PCICt(strsplit(bcca$dim$time$units, ' ')[[1]][3],
                      bcca$dim$time$calendar) + bcca$dim$time$vals*86400
dates.bcca <-  apply(do.call(rbind, strsplit(format(dates.bcca, '%Y %m %d'),
                    ' ')), 2, as.integer)

n.qpqm <- nrow(dates.qpqm)
n.bcca <- nrow(dates.bcca)

years <- unique(dates.bcca[,1])
n.chunks <- 2 ##length(years)

##----------------------------------------------------------
nlon <- 400



for (n in 1:n.chunks) {
  year <- years[n]
  print(year)
  t.st <- grep(year,dates.bcca[,1])[1]
  print(t.st)
  date.sub <- dates.bcca[dates.bcca[,1] %in% year,]
  t.cnt <- nrow(date.sub)
  mons <- date.sub[,2]

  read.time <- proc.time()
  var.bcca <- ncvar_get(bcca,varid=var.name,start=c(1,1,t.st),count=c(nlon,nlat,t.cnt))
  var.qpqm <- ncvar_get(qpqm,varid=var.name,start=c(1,1,t.st),count=c(nlon,nlat,t.cnt))
  var.final <- var.qpqm*0

  ncvar_put(nc=out, varid=var.name, vals=var.final,
            start=c(1, 1, t.st), count=dim(var.final))  
  print('Reading data time')
  print(proc.time()-read.time)
  
  print('Arrange into list')
  mn.subset <- vector(mode='list',length=12)
  bcca.list <- vector(mode='list',length=12)
  qpqm.list <- vector(mode='list',length=12)
  for (m in 1:12) {
    subset <- which(mons==m)
    mn.subset[[m]] <- subset
    bcca.list[[m]] <- var.bcca[,,subset]
    qpqm.list[[m]] <- var.qpqm[,,subset]
  }
  ##qpqm.list=qpqm.list,
  print('Parallel')
  par.time <- proc.time()
  ##for (i in 1:12) {
  dqm <- foreach(i=1:12,
                 b.list=bcca.list,
                 q.list=qpqm.list, ##, ##) %dopar% { 
                 .inorder=TRUE,
                 .options.mpi=mpi.options) %dopar% { 
                   ##.combine='cbind',
                   ##.multicombine=TRUE,
                   
                   b.list <- jitter(bcca.list[[i]],0.01)
                   q.list <- qpqm.list[[i]]
                   
                   slen <- dim(b.list)[3]
                   nlon <- dim(b.list)[1]
                   nlat <- dim(b.list)[2]
                   bcca.ranks <- round(aperm(apply(b.list,c(1,2),rank,ties.method='average'),c(2,3,1))) ##Permute to fix the rank permute  
                   data.p <- aperm(q.list,c(3,1,2))
                   
                   dim(data.p) <- c(slen,nlon*nlat)
                   data.l <- lapply(apply(data.p,2,as.list),unlist)
                   
                   index.p <- aperm(bcca.ranks,c(3,1,2))
                   dim(index.p) <- c(slen,nlon*nlat)
                   index.l <- lapply(apply(index.p,2,as.list),unlist)
                   
                   result <- mapply(reorder,data.l,index.l)
                   dim(result) <- c(slen,nlon,nlat)
                   final <- aperm(result,c(2,3,1))
                 }
  print('Done parallel')
  print(proc.time() - par.time)
  for (l in 1:12)
    var.final[,,mn.subset[[l]]] <- dqm[[l]]
  
  ncvar_put(nc=out, varid=var.name, vals=var.final,
            start=c(1, 1, t.st), count=dim(var.final))  
}


print('Elapsed time')
print(proc.time() - ptm)
#browser()  


nc_close(qpqm)
nc_close(bcca)
nc_close(out)

closeCluster(cl)
mpi.quit()

################################################################################
