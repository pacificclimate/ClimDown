##******************************************************************************
# Bias Corrected Constructed Analogue (BCCA) downscaling algorithm
# Based loosely off of code by Alex Cannon <acannon@uvic.ca>
# Rewritten by James Hiebert <hiebert@uvic.ca>

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
  obs.time <- netcdf.calendar(nc.obs, 'time')

  # Figure out which GCM grid boxes are associated with each fine-scale grid point
  grid.mapping <- regrid.coarse.to.fine(gcm.lats, gcm.lons, obs.lats, obs.lons)
  xi <- grid.mapping$xi
  yi <- grid.mapping$yi

  xn <- length(unique(as.vector(xi)))
  yn <- length(unique(as.vector(yi)))
  aggregates <- array(dim=c(length(gcm.lons), length(gcm.lats), length(obs.time)))

  chunk.size <- optimal.chunk.size(length(obs.lons) * length(obs.lats))

  chunks <- chunk.indices(length(obs.time), chunk.size)
  # Loop over chunks fo time
  for (i in chunks) {
    print(paste("Aggregating timesteps", i['start'], "-", i['stop'], "/", length(obs.time)))
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
                             historical.start=getOption('calibration.start'),
                             historical.end=getOption('calibration.end'),
                             detrend=FALSE,
                             ratio=FALSE) {
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
                             ratio=ratio, detrend=detrend, n.max=NULL)
            c(dqm.tmp$g.p.bc, dqm.tmp$g.h.bc, dqm.tmp$g.f.bc)
        }
    }, points[,'row'], points[,'col']))
    rv
}

# x: a vector of lat x lon x num.analogues
# weights: a vector of length num.analogues
apply.analogue <- function(x, weights) {
    n.cells <- prod(dim(x)[1:2])
    weights <- sapply(weights, rep, n.cells)
    dim(weights) <- dim(x)
    apply(x * weights, 1:2, sum)
}

# analog.indices: vector of time indices that correspond to the timesteps to compose together
# weights: vector of length num.analogues corresponding to the analog indices
# obs.nc: An open netcdf file containing gridded observations
apply.analogues.netcdf <- function(analog.indices, weights, obs.nc, varid='tasmax') {
    dims <- obs.nc$var[[varid]]$size[1:2]
    apply(
        array(
            mapply(function(i, w) {
                       ncvar_get(nc=obs.nc, varid=varid,
                                 start=c(1, 1, i),
                                 count=c(-1, -1, 1)) * w
                   },
                   analog.indices, weights
                   ),
            dim=dims,
            ),
        1:2, sum
    )
}

# obs.at.analogues should be a matrix (n.analogues x number of cells)
# gcm.values should a 1d vector of gcm values for each cell at the given time step
construct.analogue.weights <- function(obs.at.analogues, gcm.values) {
    n.analogue <- nrow(obs.at.analogues)
    alib <- jitter(obs.at.analogues)
    Q <- alib %*% t(alib)
    ridge <- getOption('tol') * mean(diag(Q))
    ridge <- diag(n.analogue) * ridge
    (solve(Q + ridge) %*% alib %*% as.matrix(gcm.values))[,1]
}

# times: timeseries vector of PCICt types
# today: PCICt object, a particular day of the year (year is not important) for which
# to compute a window around
# delta.days: an integer describing the size of the window on either side of today
analogue.search.space <- function(times, today,
                                  delta.days=getOption('delta.days'),
                                  t0=getOption('calibration.start'),
                                  tn=getOption('calibration.end')) {
    cal <- attr(times, 'cal')
    dpy <- if (cal == '360') 360 else 365
    jdays <- as.numeric(format(times, '%j'))
    today <- as.numeric(format(today, '%j'))

    distance <- abs(jdays - today)
    in.days <- distance <= delta.days | distance >= (dpy - delta.days)
    in.years <- times >= as.PCICt(t0, cal=cal) & times <= as.PCICt(tn, cal=cal)
    which(in.days & in.years)
}

# FIXME: There's one remaining difference between pr and tas in the old BCCA code that
# needs to be explained and implemented
# 116,11yc108,109
# <     pr.gcm.i <- pr.gcm[i,na.mask]^expon
# <     pr.agg.alib <- pr.aggregate[alib,na.mask]^expon
# ---
# >     tasmin.gcm.i <- tasmin.gcm[i,na.mask]
# >     tasmin.agg.alib <- tasmin.aggregate[alib,na.mask]
# And then when the analogues are applied...
# 133,134c125
# <         pr.analogue.j <- round(pr.analogue.j, 3)
# <         pr.analogue <- pr.analogue + pr.weights[j]*(pr.analogue.j^expon)
# ---
# >         tasmin.analogue <- tasmin.analogue + tasmin.weights[j]*(tasmin.analogue.j)
# 136,137c127
# <     pr.analogue[pr.analogue < 0] <- 0
# <     pr.analogue <- pr.analogue^(1/expon)
# ---
# >     cat('*')
# What's the purpose of taking the square root before combining the analogues and then
# squaring after they're combined?

# gcm: a 2d vector representing a single time step of GCM valuse
# agged.obs: a 3d vector (lat x lon x time) representing the aggregated observations
# times: PCICt vector of time values for the aggregated obs
# now: PCICt value of the current time step
# returns 30 indices of the timestep for the closest analog and their corresponding weights
find.analogues <- function(gcm, agged.obs, times, now, n.analogues=getOption('n.analogues')) {
    # FIXME: alib is only be defined for each julian day in any year...
    # it is 50x redundant as defined and could be precomputed
    ti <- analogue.search.space(times, now)
    agged.obs <- agged.obs[,,ti,drop=FALSE]

    # Find the n.analogue closest observations from the library
    # (obs years * (delta days * 2 + 1)) x cells
    # substract the GCM at this time step from the aggregated obs *for every library time value*
    # square that difference

    diffs <- (agged.obs - array(gcm, dim(agged.obs))) ^ 2
    diffs <- apply(diffs, 3, sum, na.rm=T)

    # Then find the 30 lowest differences
    # returns the indices for the n closest analogues
    # of this particular GCM timestep
    analogue.indices <- quickest.select(diffs, n.analogues)

    # Constructed analogue weights
    na.mask <- !is.na(agged.obs[,,1])
    obs.at.analogues <- t(matrix(agged.obs[,,analogue.indices][na.mask], ncol=n.analogues))
    weights <- construct.analogue.weights(obs.at.analogues, gcm[na.mask])

    # Remap the search space indices to full timeseries indices
    analogue.indices <- ti[analogue.indices]

    list(analogues=analogue.indices, weights=weights)
}

# gcm: a 3d vector (lat x lon x time) representing a GCM simulation
# agged.obs: a 3d vector (lat x lon x time) representing the aggregated observations
# gcm.times: PCICt vector of time values for the GCM
# obs.time: PCICt vector of time values for the aggregated obs
find.all.analogues <- function(gcm, agged.obs, gcm.times, obs.times) {
    split(
        sapply(seq_along(gcm.times), function(i) {
                   print(paste(i, '/', length(gcm.times)))
                   find.analogues(gcm[,,i], agged.obs, obs.times, gcm.times[i])
        }),
        c('indices', 'weights')
    )
}

mk.output.ncdf <- function(file.name, varname, template.nc, global.attrs=list()) {
    nc <- nc_create(file.name, template.nc$var[[varname]])
    mapply(function(name, value) {
        ncatt_put(nc, varid=0, attname=name, attval=value)
    }, names(global.attrs), global.attrs)
    nc
}

#' @title High-level NetCDF I/O wrapper for the Bias Correction Constructed Analogues (BCCA) pipeline
#'
#' @description BCCA starts by spatially aggregating high-resolution
#' gridded observations up to the scale of a GCM. Then it proceeds to
#' bias correcting the GCM based on those observations. Finally, it
#' conducts the search for temporal analogues (which is the most
#' expensive part of the operation). This involves taking each
#' timestep in the GCM and searching for the top 30 closest timesteps
#' (for some function of "close") in the gridded observations. For
#' each of the 30 closest "analogue" timesteps, BCCA records the
#' integer number of the timestep and a weight for each of the
#' analogues. These are all saved in output.file.
#' 
#' @param gcm.file Filename of GCM simulations
#' @param obs.file Filename of high-res gridded historical observations
#' @param varname Name of the NetCDF variable to downscale (e.g. 'tasmax')
#' @return A list object with two values: 'indices' and 'weights', each of which is a vector with 30 items
#'
#' @export
bcca.netcdf.wrapper <- function(gcm.file, obs.file, varname='tasmax') {
    is.pr <- varname == 'pr'

    # Read in GCM data
    nc <- nc_open(gcm.file)
    gcm <- ncvar_get(nc, varname)

    units <- ncatt_get(nc, varname, 'units')$value
    gcm <- ud.convert(gcm, units, target.units[varname])

    if (is.pr) {
        gcm[gcm < 0] < 0
    }

    gcm.time <- netcdf.calendar(nc, 'time')
    nc_close(nc)

    print("Aggregating observations to GCM scale")
    aggd.obs <- create.aggregates(obs.file, gcm.file, varname)

    if (is.pr) {
        aggd.obs[aggd.obs < 0] < 0
    }

    print("Reading time values from the aggregated observations")
    nc <- nc_open(obs.file)
    obs.time <- netcdf.calendar(nc, 'time')
    nc_close(nc)

    print("Bias correcting the aggregated observations")
    bc.gcm <- bias.correct.dqm(
        gcm, aggd.obs, obs.time, gcm.time,
        getOption('calibration.start'), getOption('calibration.end'),
        detrend=!is.pr, ratio=is.pr
    )
    print("Finding an analogous observered timestep for each GCM time step")
    find.all.analogues(bc.gcm, aggd.obs, gcm.time, obs.time)
}
