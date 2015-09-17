##******************************************************************************
# Bias Corrected Constructed Analogue (BCCA) downscaling algorithm
# Based loosely off of code by Alex Cannon <acannon@uvic.ca>
# Rewritten by James Hiebert <hiebert@uvic.ca>

library(ncdf4)
source(paste(code.dir,'DQM.R',sep=''))
source(paste(code.dir,'netcdf.calendar.R',sep=''))

# "--args gcm.file='${gcm.file}' obs.file='${obs.file}' varid='${varid}'"
args <- commandArgs(trailingOnly=TRUE)
for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
}

na.mask <- function(grid) {
    which(!is.na(grid))
}

##******************************************************************************
# Bias correct daily GCM values using detrended quantile mapping algorithm
##******************************************************************************
bias.correct.dqm <- function(gcm, aggd.obs,
                             obs.time,
                             gcm.time,
                             historical.start=as.PCICt('1951-1-1', 'gregorian'),
                             historical.end=as.PCICt('2005-12-31', 'gregorian'),
                             detrend=FALSE) {
    t0 <- historical.start
    tn <- historical.end
    prehist.period <- gcm.time < t0
    future.period <- gcm.time > tn
    hist.period.gcm <- ! (prehist.period | future.period)
    hist.period.obs <- obs.time >= t0 & obs.time <= tn

    points <- na.mask(aggd.obs[,,1])

    for(i in seq_along(points)) {
        point <- points[i]
        if (all(is.na(gcm[point,]))) {
            gcm[point,] <- rep(NA,dim(gcm)[3])
        } else {
            dqm.tmp <- mnDQM(obs.h=aggd.obs[point, hist.period],
                             gcm.h=gcm[point, hist.period],
                             gcm.f=gcm[point, future.period],
                             months.obs.h=as.numeric(format(obs.time[hist.period.obs], '%m')),
                             months.gcm.h=as.numeric(format(gcm.time[hist.period.gcm], '%m')),
                             months.gcm.f=as.numeric(format(gcm.time[future.period], '%m')),
                             gcm.p=gcm[point, prehist.period],
                             months.gcm.p=as.numeric(format(gcm.time[prehist.period], '%m')),
                             ratio=TRUE, detrend=detrend, n.max=NULL)
            gcm[,point] <- c(dqm.tmp$g.p.bc, dqm.tmp$g.h.bc, dqm.tmp$g.f.bc)
        }
    }
    gcm
}

# NetCDF I/O wrapper for bias.correct.dqm()
bias.correct.dqm.netcdf <- function(gcm.nc, obs.nc, varname='tasmax') {
    # Read in GCM data
    nc <- nc_open(gcm.nc)
    gcm <- ncvar_get(nc, varname)-273.15 # FdIXME: This is tas, but make this call udunits
    gcm.time <- netcdf.calendar(gcm.nc, 'time', pcict=TRUE)
    nc_close(nc)
    #gcm <- sweep(gcm, 2, na.mask, '*') # Pretty sure that this is unnecessary
    gcm <- round(gcm, 3) # FIXME: This depends on the variable

    # Read the aggregated obs
    nc <- nc_open(obs.nc)
    #aggd.obs <- ncvar_get(nc, varname)
    obs.time <- netcdf.calendar(obs.nc, 'time', pcict=TRUE)
    # FIXME: replace this w/ a call to create.aggregates()
    load('tasmax.aggregates.Rdata')
    aggd.obs <- aggregates

    bias.correct.dqm(gcm, aggd.obs, obs.time, gcm.time, detrend=FALSE)
}

ptm <- proc.time()
print('Starting')

gcm <- bias.correct.dqm.netcdf(gcm.file, obs.file, varname)
save(gcm, file=paste('gcm_bc.Rdata'))

print('Elapsed time')
print(proc.time()-ptm)
