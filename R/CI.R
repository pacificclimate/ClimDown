##******************************************************************************
# Climate Imprint (CI) downscaling algorithm
# Conceptually based (very) loosely off of code by Alex Cannon <acannon@uvic.ca>
# Rewritten by James Hiebert <hiebert@uvic.ca>

is.range.subset <- function(inner, outer) {
    if (any(is.na(c(inner, outer)))) {
        FALSE
    } else {
        outer[1] < inner[1] & outer[2] > inner[2]
    }
}

# O(n) time, O(n) space
monthly.climatologies <- function(gcm, gcm.times, varname) {
    if (varname == 'pr') { # PRISM climatologies given as monthly total climatologies
        # Calculate total precipitation by month for full baseline series
        monthly.ts.factor <- factor(format(gcm.times, '%Y-%m'))
        clim.mon.factor <- factor(format(as.Date(paste0(levels(monthly.ts.factor), '-01')), '%m'))
        monthly.totals <- apply(gcm, 1:2, function(x, fac) {
            tapply(x, fac, sum, na.rm=T)
        }, monthly.ts.factor)
        # Take average of monthly totals by calendar month
        rv <- apply(monthly.totals, 2:3, function(x, fac) {
            tapply(x, fac, mean, na.rm=T)
        }, clim.mon.factor)
    } else {
        monthly.factor <- factor(format(gcm.times, '%m'))
        rv <- apply(gcm, 1:2, function(x) {
            tapply(x, monthly.factor, mean)
        })
    }
    aperm(rv, c(2, 3, 1))
}

# O(2n) time, O(2n) space
daily.anomalies <- function(gcm, gcm.times, cal.start, cal.end, varname) {
    `%op%` <- ifelse (varname == 'pr', `/`, `-`)
    ti <- compute.time.overlap(gcm.times, cal.start, cal.end)
    clima <- monthly.climatologies(gcm[,,ti], gcm.times[ti], varname)
    months <- as.integer(format(gcm.times, '%m'))
    array(
        mapply(
            function(ti, m) {
                gcm[,,ti] %op% clima[,,m]
            },
            seq_along(gcm.times), months
        ),
        dim=dim(gcm)
    )
}

# x is a 3d array; clima is a 3d array with 3rd dim of length 12;
# monthly factor is of length num_timesteps
# O(n) time, O(n) space
apply.climatologies <- function(x, clima, monthly.factor) {
    array(
        sapply(
            seq_along(monthly.factor),
            function(i) {
                t <- monthly.factor[i]
                x[,,i] + clima[,,t]
            }
        ),
        dim=dim(x)
    )
}

apply.climatologies.nc <- function(nc, clima, monthly.factor, varname='tasmax', nt.per.chunk=100) {
    `%op%` <- ifelse(varname == 'pr', `*`, `+`)
    nt <- nc$var[[varname]]$varsize[3]
    chunks <- chunk.indices(nt, nt.per.chunk)
    for (i in chunks) {
        print(paste("Applying climatologies to file", nc$filename, "steps", i['start'], ':', i['stop'], '/', nt))
        x <- CD_ncvar_get(nc, varid=varname, start=c(1, 1, i['start']), count=c(-1, -1, i['length']))
        t <- monthly.factor[i['start']:i['stop']]
        x <- x %op% clima[,,t]
        ncvar_put(nc, varid=varname, vals=x, start=c(1, 1, i['start']), count=c(-1, -1, i['length']))
        rm(x)
        gc()
    }
    NULL
}

# Computes a mean across time (3rd dimension)
# Return value is a grid the size of the first 2 dimensions of nc['var.id']
# with the 3rd dimension being the number of factors
# O(n) time, O(x*y*12 + x*y*nt.per.chunk) space
chunked.factored.running.mean <- function(nc, fact, var.id, nt.per.chunk=100) {
    dims <- nc$var[[var.id]]$varsize
    nt <- dims[3]
    # calculate how many timesteps we can get at one time

    # initialize return array and count of time steps
    rv <- array(0, dim=c(dims[1:2], length(levels(fact))))
    t <- array(0, length(levels(fact)))
    chunks <- chunk.indices(nt, nt.per.chunk)

    for (i in chunks) {
        f <- fact[i['start']:i['stop']]

        # fetch the data
        print(paste("Reading timesteps", i['start'], ':', i['stop'], '/', nt, 'from file:', nc$filename))
        x <- CD_ncvar_get(nc, varid=var.id, start=c(1, 1, i['start']), # get obs for one chunk
                       count=c(-1, -1, i['length']))

        print("Computing the temporal mean")
        # get the sum across time (to be averaged)
        subsum <- aperm(
            apply(x, 1:2, function(y) {
                tapply(y, f, sum, na.rm=FALSE)
            }), c(2, 3, 1)
        )

        lengths <- tapply(f, f, length)
        # Only work on means for months where we have values in this chunk
        l <- which(!is.na(lengths))
        lengths[is.na(lengths)] <- 0

        # increment the number of steps
        t <- t + lengths

        # compute a running mean based on the number of steps traversed so far
        previous.mean <- array(
            apply(rv[,,l,drop=F], 1:2, '*', ((t[l] - lengths[l]) / t[l])),
            dim=c(length(l), dims[1:2])
        )
        current.mean <- array(
            apply(subsum[,,l,drop=F], 1:2, '/', t[l]),
            dim=c(length(l), dims[1:2])
        )            
        rv[,,l] <- aperm( previous.mean + current.mean, c(2, 3, 1))
    }
    rv
}

xy.grid <- function(lats, lons) {
    yy <- matrix(lats, ncol=length(lats), nrow=length(lons), byrow=T)
    xx <- matrix(lons, ncol=length(lats), nrow=length(lons), byrow=F)
    list(x=xx, y=yy)
}

# O(obs_x * obs_y * gcm_t) time, O(gcm + obs_x * obs_y * gcm_t) space
interpolate.gcm.to.obs <- function(gcm.lats, gcm.lons, obs.lats, obs.lons, gcm) {
    nt <- dim(gcm)[3]
    obs.grid <- xy.grid(obs.lats, obs.lons)
    src <- list(x=gcm.lons, y=gcm.lats)
    dst <- matrix(c(obs.grid$x, obs.grid$y), ncol=2)
    array(
         apply(gcm, 3, function(z) {
            cat('*')
            src$z <- z
            interp.surface(src, dst) # from fields pacakge
        }), dim=c(dim(obs.grid$x), nt)
    )
}

chunked.interpolate.gcm.to.obs <- function(gcm.lats, gcm.lons,
                                           obs.lats, obs.lons,
                                           gcm, output.nc, varname,
                                           nt.per.chunk=100) {
    nt <- dim(gcm)[3]
    ncells <- length(gcm.lats) * length(gcm.lons)
    print(output.nc$var[[varname]]$varsize)
    print(c(length(obs.lons), length(obs.lats), dim(gcm)[3]))
    stopifnot(output.nc$var[[varname]]$varsize == c(length(obs.lons), length(obs.lats), dim(gcm)[3]))
    obs.grid <- xy.grid(obs.lats, obs.lons)

    if (!(is.range.subset(range(obs.lons), range(gcm.lons)) &
          is.range.subset(range(obs.lats), range(gcm.lats)))) {
        stop("Observation domain must be a proper spatial subset of the GCM domain (but it's not). Please check your input files.")
    }

    src <- list(x=gcm.lons, y=gcm.lats)
    dst <- matrix(c(obs.grid$x, obs.grid$y), ncol=2)

    chunks <- chunk.indices(nt, nt.per.chunk)

    for (i in chunks) {
        i0 <- i['start']
        iN <- i['stop']
        print(paste("Interpolating timesteps", i0, "-", iN, "/", nt, "to file", output.nc$filename))
        rv <- array(
            apply(gcm[,,i0:iN], 3, function(z) {
                src$z <- z
                interp.surface(src, dst)
            }),
            dim=c(length(obs.lons), length(obs.lats), i['length'])
        )
        ncvar_put(output.nc, varname, vals=rv, start=c(1, 1, i0), count=c(-1, -1, i['length']))
        rm(rv)
        gc()
    }
    nc_sync(output.nc)
    NULL
}

# FIXME: this name is duplicated from CA.R
mk.output.ncdf <- function(file.name, gcm.varname, obs.varname, gcm.template, obs.template, global.attrs=list()) {
    dims <- c(obs.template$var[[obs.varname]]$dim[1:2], gcm.template$var[[gcm.varname]]$dim[3])
    var <- ncvar_def(gcm.varname, getOption('target.units')[gcm.varname], dims, NA)
    nc <- nc_create(file.name, var)
    mapply(function(name, value) {
        ncatt_put(nc, varid=0, attname=name, attval=value)
    }, names(global.attrs), global.attrs)
    nc
}

nc_getx <- function(nc) {
    lons <- ncvar_get(nc, 'lon')
    if (max(lons) > 180) {
        lons <- lons - 360
    }
    lons
}

nc_gety <- function(nc) {
    sort(ncvar_get(nc, 'lat'))
}

is.clim.file <- function(obs) {
    "climatology_bnds" %in% names(obs$var) | "climatology_bounds" %in% names(obs$var)
}

#' @title High-level NetCDF wrapper for Climate Imprint (CI)
#'
#' @description CI performs several steps. For the GCM input it
#' calculates daily climate anomalies from a given calibration period
#' (default 1951-2005). These daily GCM anomalies are interpolated to
#' the high-resolution observational grid. These interpolated daily
#' anomalies constitute the "Climate Imprint". If the high resolution
#' gridded observations are daily observations instead of pre-computed
#' monthly climatologies, they are then grouped into months and a climatology
#' is calculated for each month. Finally the observed climatology is
#' added to the GCM-based climate imprint and the final result is
#' saved to output.file.
#'
#' @param gcm.file Filename of GCM simulations
#' @param obs.file Filename of high-res gridded historical observations 
#' (either daily data or monthly climatologies)
#' @param output.file Filename to create (or overwrite) with the climate imprint outputs
#' @param varname Name of the NetCDF variable to downscale (e.g. 'tasmax')
#'
#' @examples
#' \dontrun{
#' ci.nc <- tempfile(fileext='.nc')
#' ClimDown::ci.netcdf.wrapper('./tiny_gcm.nc', './tiny_obs.nc', ci.nc)
#' unlink(ci.nc)
#' }
#'
#' @references Hunter, R. D., & Meentemeyer, R. K. (2005). Climatologically aided mapping of daily precipitation and temperature. Journal of Applied Meteorology, 44(10), 1501-1510.
#'
#' Ahmed, K. F., Wang, G., Silander, J., Wilson, A. M., Allen, J. M., Horton, R., & Anyah, R. (2013). Statistical downscaling and bias correction of climate model outputs for climate change impact assessment in the US northeast. Global and Planetary Change, 100, 320-332.
#' @export
ci.netcdf.wrapper <- function(gcm.file, obs.file, output.file) {

    gcm.varname <- getOption('gcm.varname')
    obs.varname <- getOption('obs.varname')
    
    nc.gcm <- nc_open(gcm.file)
    gcm <- CD_ncvar_get(nc.gcm, gcm.varname)
    gcm.lats <- nc_gety(nc.gcm)
    gcm.lons <- nc_getx(nc.gcm)
    gcm.times <- netcdf.calendar(nc.gcm)

    print('Calculating daily anomalies on the GCM')
    cal <- attr(gcm.times, 'cal')
    cstart <- as.PCICt(getOption('calibration.start'), cal=cal)
    cend <- as.PCICt(getOption('calibration.end'), cal=cal)
    anom <- daily.anomalies(gcm, gcm.times, cstart, cend, gcm.varname)

    nc.obs <- nc_open(obs.file)
    obs.lats <- nc_gety(nc.obs)
    obs.lons <- nc_getx(nc.obs)
    obs.times <- netcdf.calendar(nc.obs)

    # Calculate the chunk size
    nt.per.chunk <- optimal.chunk.size(length(obs.lats) * length(obs.lons), getOption('max.GB'))

    print('Creating cache file for the interpolated GCM')
    output.nc <-mk.output.ncdf(output.file, gcm.varname, obs.varname, nc.gcm, nc.obs)
    nc_close(nc.gcm)

    print('Interpolating the GCM daily anomalies to observation grid')
    chunked.interpolate.gcm.to.obs(gcm.lats, gcm.lons, obs.lats, obs.lons, anom, output.nc, gcm.varname, nt.per.chunk)

    print('Check observations file')
    if (is.clim.file(nc.obs)) {
        print('Reading the monthly climatologies from the observations')
        monthly.climatologies <- ncvar_get(nc.obs, obs.varname)
    } else {
        print('Calculating the monthly factor across the observation time series')
        monthly.factor <- factor(as.numeric(format(obs.times, '%m')))
        
        print('Calculating the monthly climatologies for the observations')
        monthly.climatologies <- chunked.factored.running.mean(nc.obs, monthly.factor, obs.varname, nt.per.chunk)   
    }
    
    nc_close(nc.obs)

    print('Calculating the monthly factor across the GCM time series')
    monthly.factor <- factor(as.numeric(format(gcm.times, '%m')))

    print('Adding the monthly climatologies to the interpolated GCM')
    apply.climatologies.nc(output.nc, monthly.climatologies, monthly.factor, gcm.varname, nt.per.chunk)
    nc_close(output.nc)
}

