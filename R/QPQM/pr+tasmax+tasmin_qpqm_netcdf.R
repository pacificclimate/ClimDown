#!/usr/bin/Rscript

################################################################################

ptm <- proc.time()

args <- commandArgs(TRUE)

#args <- c('/home/data/climate/downscale/CMIP5/nobackup/QPQM/nobackup/pr+tasmax+tasmin_day_ANUSPLIN300_observation_v20140414_19710101-20001231.nc',
#          '/home/data/climate/downscale/CMIP5/nobackup/QPQM/done/pr+tasmax+tasmin_day_QPQM+ANUSPLIN300+ACCESS1-0_historical+rcp45_r1i1p1_19500101-21001231.nc',
#          '/datasets/projects-rci-scratch/test/dqm_netcdf.nc','2')

## Help section
if('--help' %in% args || length(args)==0) {
    cat('
      Quantile perturbation quantile mapping

      Usage:
      ./pr+tasmax+tasmin_qpqm_netcdf.R [obs_file] [gcm_file] [output_file] [cores]

      Arguments:
      obs_file - Gridded historical observations
      gcm_file - GCM simulations interpolated to the obs_file grid
      output_file - The file to create (or overwrite) with the bias corrected outputs
      varname - The name of the variable to downscale
      cores - The number of processor cores

      All files should have the same spatial domain.\n\n')

      q(save='no')
}

library(ncdf4)
library(abind)
library(PCICt)
#library(Rmpi)
#library(doMPI)	
source('/home/hiebert/code/git/ClimDown/R/QPQM/QPQM.R')

#cl <- startMPIcluster(count=args[4])
#registerDoMPI(cl)
#exportDoMPI(cl, c('QPQM', 'tQPQM'))
#mpi.options <- list(info=FALSE, profile=FALSE)
#cat('cores =', getDoParWorkers(), '\n')

n.window <- 1
multiyear <- TRUE
expand.multiyear <- TRUE
n.multiyear <- 30
trace <- 0.005
jitter.factor <- 0.01
cstart <- 1971 
cend <- 2000
n.chunks <- 500

pr.n.tau <- 1001
tasmax.n.tau <- 101
tasmin.n.tau <- 101

seasonal <- list(pr=TRUE, tasmax=TRUE, tasmin=TRUE)
ratio <- seasonal

pr.offset <- 0.
pr.scale <- 1.
tasmax.offset <- 0.
tasmax.scale <- 1.
tasmin.offset <- 0.
tasmin.scale <- 1.

# Read in the input and output files
obs.file <- args[1]
gcm.file <- args[2]
out.file <- args[3]
varname <- args[4]

gcm <- nc_open(gcm.file)
obs <- nc_open(obs.file)

cat('Copying input file', gcm.file, 'to', out.file, '\n')
#file.copy(from=gcm.file, to=out.file, overwrite=TRUE)
out <- nc_create(out.file, obs$var[[varname]])

lat <- gcm$dim$lat$vals
lon <- gcm$dim$lon$vals

dates.gcm <- as.PCICt(strsplit(gcm$dim$time$units, ' ')[[1]][3],
                      gcm$dim$time$calendar) + gcm$dim$time$vals*86400
dates.gcm <-  apply(do.call(rbind, strsplit(format(dates.gcm, '%Y %m %d'),
                    ' ')), 2, as.integer)

dates.obs <- as.PCICt(strsplit(obs$dim$time$units, ' ')[[1]][3],
                      obs$dim$time$calendar) + obs$dim$time$vals*86400
dates.obs <-  apply(do.call(rbind, strsplit(format(dates.obs, '%Y %m %d'),
                    ' ')), 2, as.integer)

cal.gcm <- c(min(which(dates.gcm[,1]==cstart)), max(which(dates.gcm[,1]==cend)))
cal.obs <- c(min(which(dates.obs[,1]==cstart)), max(which(dates.obs[,1]==cend)))

n.gcm <- nrow(dates.gcm)
n.obs <- diff(cal.obs)+1

dates.o.c <- dates.obs[cal.obs[1]:(cal.obs[1]+diff(cal.obs)),]
dates.m.c <- dates.gcm[cal.gcm[1]:(cal.gcm[1]+diff(cal.gcm)),]
na.gcm <- rep(NA, n.gcm)

# A very strange way of dividing up into equal chunks
i.indices <- seq_along(lon)
i.chunks <- suppressWarnings(matrix(c(i.indices, rep(NA, length(i.indices))),
                             ncol=n.chunks*2))
i.chunks <- i.chunks[,1:(min(which(is.na(i.chunks[1,])))-1), drop=FALSE]
n.chunks <- ncol(i.chunks)

jj <- seq_along(lat)

browser()

for(chunk in 1:2) { ##ncol(i.chunks)){
    ii <- na.omit(i.chunks[,chunk])
    ## I'm pretty sure that this chunking is backwards from optimal. The larger the ii the better
    ## Experiment measure 1068 x 10 vs 10 x 1068 (each are just under a GB)
    ## Actualy, they seem to be about the same, so, whatever
    cat('--> bias correcting', varname, 'chunk', chunk, '/', n.chunks, '-')
    o.c.chunk <- ncvar_get(obs, start=c(ii[1], jj[1], cal.obs[1]),
                           count=c(length(ii), length(jj), n.obs),
                           varid=varname, collapse_degen=FALSE)
    cat('-')
    m.p.chunk <- ncvar_get(gcm, start=c(ii[1], jj[1], 1),
                           count=c(length(ii), length(jj), n.gcm),
                           varid=varname, collapse_degen=FALSE)

    # FIXME: dont bother permuting and redimming... just do an apply(preserve=3) below
    m.p.chunk <- aperm(m.p.chunk, c(3, 1, 2))
    dim(m.p.chunk) <- c(n.gcm, length(ii) * length(jj))
    o.c.chunk <- aperm(o.c.chunk, c(3, 1, 2))
    dim(o.c.chunk) <- c(n.obs, length(ii) * length(jj))
    cat('*\n')
    ij <- t(expand.grid(i=seq_along(ii), j=seq_along(jj)))

    m.p.chunk <- foreach(ij=ij,
                         o.c=apply(o.c.chunk, 2, identity),
                         m.p=apply(m.p.chunk, 2, identity),
                         .inorder=TRUE,
                         .combine=cbind,
                         .multicombine=TRUE
                         #.options.mpi=mpi.options) %do% {
                         ) %do {
        i = ij['i',]
        j = ij['j',]
        if(all(is.na(o.c)) || all(is.na(m.p))) {
            na.gcm
        } else {
            m.c <- m.p[cal.gcm[1]:(cal.gcm[1]+diff(cal.gcm))]
            tQPQM(o.c=o.c, m.c=m.c, m.p=m.p, dates.o.c=dates.o.c,
                  dates.m.c=dates.m.c, dates.m.p=dates.gcm,
                  n.window=n.window, ratio=ratio[[varname]], trace=trace,
                  jitter.factor=jitter.factor, seasonal=seasonal[[varname]],
                  multiyear=multiyear, n.multiyear=n.multiyear,
                  expand.multiyear=expand.multiyear, n.tau=pr.n.tau)
        }
    }
    dim(m.p.chunk) <- c(n.gcm, length(ii), length(jj))
    cat(dim(m.p.chunk))
    m.p.chunk <- aperm(m.p.chunk, c(2, 3, 1))
    cat('--> writing chunk', chunk, '/', n.chunks, '\n')
    ncvar_put(nc=out, varid='pr', vals=m.p.chunk,
              start=c(ii[1], jj[1], 1), count=dim(m.p.chunk))
    print(object.size(x=lapply(ls(), get)), units="Mb")
    gc()
}

nc_close(gcm)
nc_close(obs)
nc_close(out)

print('Elapsed Time')
print(proc.time() - ptm)

#closeCluster(cl)
#mpi.quit()

################################################################################
