##******************************************************************************
# Bias Corrected Constructed Analogue (BCCA) downscaling algorithm
# Based loosely off of code by Alex Cannon <acannon@uvic.ca>
# Rewritten by James Hiebert <hiebert@uvic.ca>

library(ncdf4)
code.dir <- Sys.getenv('CODE_DIR')
source(paste(code.dir, 'bisect.R', sep='/'))
source(paste(code.dir, 'netcdf.calendar.R', sep='/'))
source(paste(code.dir,'DQM.R',sep='/'))

options(max.GB=1)
options(trimmed.mean=0)
options(delta.days=45)

# "--args gcm.file='${gcm.file}' obs.file='${obs.file}' varid='${varid}'"
args <- commandArgs(trailingOnly=TRUE)
for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
}

# Input is a factor of which maps fine-scale obs cells to large-scale gcm cells
# The factor should be of length x * y
# and a 3d array of obs (x, y, time)
aggregate.obs <- function(cell.factor, obs) {
  trim <- getOption('trimmed.mean')
  apply(obs, 3, function(x) tapply(x, cell.factor, mean, trim=trim, na.rm=TRUE))
}

# Input cell indicies mapping obs grid to GCM grid
# and a 3d array of obs (x, y, time)
# Output is 3d array, gcmx x gcm y x time
aggregate.obs.to.gcm.grid <- function(xi, yi, xn, yn, obs) {
  cell.number <- xi * max(yi) + yi
  cell.factor <- factor(cell.number, unique(as.vector(cell.number)))
  # apply preserving time (dim 3)
  # so for each time step, aggregate (take the mean) according to the cell map
  rv <- aggregate.obs(cell.factor, obs)
  ti <- dim(obs)[3]
  dim(rv) <- c(xn, yn, ti)
  return(rv)
}

# Takes a vector length and chunk size
# returns a list of (start, stop, length)
chunk.indices <- function(total.size, chunk.size) {
  lapply(
    split(1:total.size, ceiling(1:total.size / chunk.size)),
    function(x) {c('start'=min(x), 'stop'=max(x), 'length'=length(x))}
    )
}

optimal.chunk.size <- function(n.elements, max.GB=getOption('max.GB')) {
  # 8 byte numerics
  floor(max.GB * 2 ** 30 / 8 / n.elements)
}

##******************************************************************************
# Read fine-scale grid and spatially aggregate to GCM grid
##******************************************************************************
create.aggregates <- function(obs.file, gcm.file, varid) {

  # Read fine-scale and GCM grid dimensions
  nc.obs <- nc_open(obs.file)
  nc.gcm <- nc_open(gcm.file)
  obs.lons <- ncvar_get(nc.obs, 'lon')
  obs.lats <- ncvar_get(nc.obs, 'lat')
  gcm.lons <- ncvar_get(nc.gcm, 'lon')-360
  gcm.lats <- ncvar_get(nc.gcm, 'lat')
  obs.time <- netcdf.calendar(nc.obs)

  # Figure out which GCM grid boxes are associated with each fine-scale grid point
  grid.mapping <- regrid.coarse.to.fine(gcm.lats, gcm.lons, obs.lats, obs.lons)
  xi <- grid.mapping$xi
  yi <- grid.mapping$yi

  xn <- length(unique(as.vector(xi)))
  yn <- length(unique(as.vector(yi)))
  aggregates <- array(dim=c(length(gcm.lons), length(gcm.lats), length(obs.time)))

  chunk.size <- optimal.chunk.size(length(obs.lons) * length(obs.lats))

  chunks <- chunk.indices(nrow(obs.time), chunk.size)
  # Loop over chunks fo time
  for (i in chunks) {
    cat(i['start'], i['stop'], '\n')
    obs <- ncvar_get(nc.obs, varid=varid, start=c(1, 1, i['start']), # get obs for one chunk
                     count=c(-1, -1, i['length']))
    agg <- aggregate.obs.to.gcm.grid(xi, yi, xn, yn, obs)
    aggregates[min(xi):max(xi), min(yi):max(yi), i['start']:i['stop']] <- agg
    rm(obs)
    gc()
  }
  aggregates
}

na.unmasked <- function(grid) {
    which(!is.na(grid), arr.ind=TRUE)
}

na.masked <- function(grid) {
    which(is.na(grid), arr.ind=TRUE)
}

##******************************************************************************
# Bias correct daily GCM values using detrended quantile mapping algorithm
##******************************************************************************
bias.correct.dqm <- function(gcm, aggd.obs,
                             obs.time,
                             gcm.time,
                             historical.start='1951-1-1',
                             historical.end='2005-12-31',
                             detrend=FALSE) {
    t0 <- as.PCICt(historical.start, attr(gcm.time, 'cal'))
    tn <- as.PCICt(historical.end, attr(gcm.time, 'cal'))
    prehist.period <- gcm.time < t0
    future.period <- gcm.time > tn
    hist.period.gcm <- ! (prehist.period | future.period)

    t0 <- as.PCICt(historical.start, attr(obs.time, 'cal'))
    tn <- as.PCICt(historical.end, attr(obs.time, 'cal'))
    hist.period.obs <- obs.time >= t0 & obs.time <= tn

    points <- na.unmasked(aggd.obs[,,1])

    rv <- array(dim=dim(gcm))
    rv[!is.na(aggd.obs[,,1])] <- t(
    mapply(function(x, y) {
        if (all(is.na(gcm[x,y,]))) {
            rep(NA,dim(gcm)[3])
        } else {
            dqm.tmp <- mnDQM(obs.h=aggd.obs[x,y, hist.period.obs],
                             gcm.h=gcm[x,y, hist.period.gcm],
                             gcm.f=gcm[x,y, future.period],
                             months.obs.h=as.numeric(format(obs.time[hist.period.obs], '%m')),
                             months.gcm.h=as.numeric(format(gcm.time[hist.period.gcm], '%m')),
                             months.gcm.f=as.numeric(format(gcm.time[future.period], '%m')),
                             gcm.p=gcm[x,y, prehist.period],
                             months.gcm.p=as.numeric(format(gcm.time[prehist.period], '%m')),
                             ratio=FALSE, detrend=detrend, n.max=NULL) # FIXME: ratio and detrend are based on variable
            c(dqm.tmp$g.p.bc, dqm.tmp$g.h.bc, dqm.tmp$g.f.bc)
        }
    }, points[,'row'], points[,'col']))
    rv
}

# NetCDF I/O wrapper for bias.correct.dqm()
bias.correct.dqm.netcdf <- function(gcm.nc, obs.nc, varname='tasmax') {
    # Read in GCM data
    nc <- nc_open(gcm.nc)
    gcm <- ncvar_get(nc, varname)-273.15 # FIXME: This is tas, but make this call udunits
    gcm.time <- netcdf.calendar(nc, 'time', pcict=TRUE)
    nc_close(nc)
    #gcm <- sweep(gcm, 2, na.mask, '*') # Pretty sure that this is unnecessary
    gcm <- round(gcm, 3) # FIXME: This depends on the variable

    # Read the aggregated obs
    nc <- nc_open(obs.nc)
    obs.time <- netcdf.calendar(nc, 'time', pcict=TRUE)
    nc_close(nc)

    print(paste('Starting spatial aggregation', obs.file, gcm.file, varid))
    ptm <- proc.time()
    aggd.obs <- create.aggregates(obs.nc, gcm.nc, varname)
    print('Elapsed time')
    print(proc.time() - ptm)

    print(paste('Starting bias correction', gcm.nc, obs.nc, varname))
    ptm <- proc.time()
    gcm <- bias.correct.dqm(gcm, aggd.obs, obs.time, gcm.time, detrend=FALSE)
    print('Elapsed time')
    print(proc.time()-ptm)
    gcm
}

# x: a vector of lat x lon x num.analogues
# weights: a vector of length num.analogues
apply.analogue <- function(x, weights) {
    n.cells <- prod(dim(x)[1:2])
    weights <- sapply(weigths, rep, n.cells)
    dim(weights) <- dim(x)
    apply(x * weights, 1:2, sum)
}

# analog.indices: vector of time indices that correspond to the timesteps to compose together
# weights: vector of length num.analogues corresponding to the analog indices
# obs.nc: An open netcdf file containing gridded observations
apply.analogues.netcdf <- function(analog.indices, weights, obs.nc) {
    sum(
        mapply(function(i, w) {
            ncvar_get(nc=obs.nc, varid=varid,
                      start=c(1, 1, i)
                      count=c(-1, -1, 1)) * w
        }
        )
    )
}

# obs.at.analogues should be a matrix (n.analogues x number of cells)
# gcm.values should a 1d vector of gcm values for each cell at the given time step
construct.analogue.weights <- function(obs.at.analogues, gcm.values) {
    n.analogue <- nrow(obs.at.analogues)
    alib <- jitter(obs.at.analogues)
    Q <- alib %*% t(alib)
    ridge <- tol * mean(diag(Q))
    ridge <- diag(n.analogue) * ridge
    solve(Q + ridge) %*% alib %*% as.matrix(gcm.values)
}

# times: timeseries vector of PCICt types
# today: PCICt object, a particular day of the year (year is not important) for which
# to compute a window around
# delta.days: an integer describing the size of the window on either side of today
analogue.search.space <- function(times, today,
                                  delta.days=getOption('delta.days'),
                                  year.range=c(1951, 2005)) {
    dpy <- attr(times, 'dpy')
    cal <- attr(times, 'cal')
    jdays <- as.numeric(strftime(times, '%j'))
    today <- as.numeric(strftime(today, '%j'))

    distance <- abs(jdays %% dpy - today)
    in.days <- distance <= delta.days | distance >= (dpy - delta.days)
    in.years <- times >= as.PCICt(paste(year.range[1], 1, 1, sep='-'), cal=cal) &
        times <= as.PCICt(paste(year.range[2], 12, 31, sep='-'), cal=cal)
    which(in.days & in.years)
}

# gcm: a 2d vector representing a single time step of GCM valuse
# agged.obs: a 3d vector (lat x lon x time) representing the aggregated observations
# times: PCICt vector of time values for the aggregated obs
# now: PCICt value of the current time step
# returns 30 indices of the timestep for the closest analog and their corresponding weights
find.analogues <- function(gcm, agged.obs, times, now) {
    # FIXME: alib is only be defined for each julian day in any year...
    # it is 50x redundant as defined and could be precomputed
    ti <- analogue.search.space(times, today)
    agged.obs <- aggd.obs[,,ti,drop=FALSE]

    # Find the n.analogue closest observations from the library
    # (obs years * (delta days * 2 + 1)) x cells
    # substract the GCM at this time step from the aggregated obs *for every library time value*
    # square that difference and then find the 30 lowest differences
    # returns n analogues for this particular GCM timestep
    diffs <- sweep(agged.obs, 1:2, gcm, '-')^2

    # FIXME, replace which(rank) with sort()[1:30]
    analogue.indices <- which(rank(apply(diffs, 3, sum, na.rm=T), ties.method='random') <= n.analogues)

    # Constructed analogue weights
    na.mask <- !is.na(aggd.obs[,,1])
    obs.at.analogues <- t(matrix(agged.obs[,,analogues][na.mask], ncol=n.analogues))
    weights <- construct.analogue.weights(obs.at.analogues, gcm[na.mask])

    list(analogues=analogue.indices, weights=weigths)
}


gcm <- bias.correct.dqm.netcdf(gcm.file, obs.file, varid)
save(gcm, file=paste('gcm_bc.Rdata'))
