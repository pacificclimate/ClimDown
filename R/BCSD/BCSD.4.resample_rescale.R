##******************************************************************************
# Bias Corrected Spatial Disaggregation (BCSD) downscaling algorithm
# Alex Cannon (acannon@uvic.ca)
##******************************************************************************
# Spatial disaggregation of GCM data, monthly resampling, and day-by-day
# scaling to match monthly means
##******************************************************************************

rm(list=ls())
library(ncdf4)
library(fields)

config <- commandArgs(trailingOnly=TRUE)
if(length(config)==0) config <- 'BCSD.config'
print(readLines(config))
source(config)

nc.obs <- nc_open(nc.obs.file)

system(paste('cp ', template.file, ' ', output.dir, output.file,
       output.suffix, '.nc', sep=''))
nc.bcsd <- nc_open(paste(output.dir, output.file, output.suffix, '.nc',
                   sep=''), write=TRUE)

ncatt_put(nc.bcsd, varid=0, attname='title', attval=output.nc.title)
ncatt_put(nc.bcsd, varid=0, attname='institution',
          attval=output.nc.institution)
ncatt_put(nc.bcsd, varid=0, attname='source', attval=output.nc.source)
ncatt_put(nc.bcsd, varid=0, attname='input_data', attval=output.nc.input_data)
ncatt_put(nc.bcsd, varid=0, attname='reference', attval=output.nc.reference)
ncatt_put(nc.bcsd, varid=0, attname='project_id', attval=output.nc.project_id)
ncatt_put(nc.bcsd, varid=0, attname='experiment_id',
          attval=output.nc.experiment_id)
ncatt_put(nc.bcsd, varid=0, attname='version', attval=output.nc.version)
ncatt_put(nc.bcsd, varid=0, attname='version_comment',
          attval=output.nc.version_comment)
ncatt_put(nc.bcsd, varid=0, attname='contact1', attval=output.nc.contact1)
ncatt_put(nc.bcsd, varid=0, attname='contact2', attval=output.nc.contact2)
ncatt_put(nc.bcsd, varid=0, attname='contact3', attval=output.nc.contact3)
ncatt_put(nc.bcsd, varid=0, attname='history', attval=output.nc.history)
ncatt_put(nc.bcsd, varid='time', attname='units', attval=output.nc.time.units)
ncatt_put(nc.bcsd, varid='time', attname='calendar', attval='gregorian')
nc_sync(nc.bcsd)

load(paste(output.dir, 'obs.time', output.suffix, '.RData', sep=''))
load(paste(output.dir, 'obs.lons', output.suffix, '.RData', sep=''))
load(paste(output.dir, 'obs.lats', output.suffix, '.RData', sep=''))

load(paste(output.dir, 'pr.obs.clim', output.suffix, '.RData', sep=''))
load(paste(output.dir, 'pr.obs.max', output.suffix, '.RData', sep=''))

load(paste(output.dir, 'tasmax.obs.clim', output.suffix, '.RData', sep=''))
load(paste(output.dir, 'tasmin.obs.clim', output.suffix, '.RData', sep=''))
load(paste(output.dir, 'wind.obs.clim', output.suffix, '.RData', sep=''))

load(paste(output.dir, 'gcm.lats', output.suffix, '.RData', sep=''))
load(paste(output.dir, 'gcm.lons', output.suffix, '.RData', sep=''))

load(paste(output.dir, 'obs.time.mn', output.suffix, '.RData', sep=''))
load(paste(output.dir, 'pr.aggregate.mn', output.suffix, '.RData', sep=''))
load(paste(output.dir, 'tasmax.aggregate.mn', output.suffix, '.RData', sep=''))
load(paste(output.dir, 'tasmin.aggregate.mn', output.suffix, '.RData', sep=''))
load(paste(output.dir, 'wind.aggregate.mn', output.suffix, '.RData', sep=''))

load(paste(output.dir, 'gcm.time.mn', output.suffix, '.RData', sep=''))
load(paste(output.dir, 'pr.gcm_bc.mn', output.suffix, '.RData', sep=''))
load(paste(output.dir, 'tasmax.gcm_bc.mn', output.suffix, '.RData', sep=''))
load(paste(output.dir, 'tasmin.gcm_bc.mn', output.suffix, '.RData', sep=''))
load(paste(output.dir, 'wind.gcm_bc.mn', output.suffix, '.RData', sep=''))

################################################################################

pr.missing_value <- ncatt_get(nc.bcsd, varid='pr',
                              attname='missing_value')$value
tasmax.missing_value <- ncatt_get(nc.bcsd, varid='tasmax',
                                  attname='missing_value')$value
tasmin.missing_value <- ncatt_get(nc.bcsd, varid='tasmin',
                                  attname='missing_value')$value
wind.missing_value <- ncatt_get(nc.bcsd, varid='wind',
                                attname='missing_value')$value

pr.scale_factor <- ncatt_get(nc.bcsd, varid='pr',
                             attname='scale_factor')$value
tasmax.scale_factor <- ncatt_get(nc.bcsd, varid='tasmax',
                                 attname='scale_factor')$value
tasmin.scale_factor <- ncatt_get(nc.bcsd, varid='tasmin',
                                 attname='scale_factor')$value
wind.scale_factor <- ncatt_get(nc.bcsd, varid='wind',
                               attname='scale_factor')$value

pr.add_offset <- ncatt_get(nc.bcsd, varid='pr',
                           attname='add_offset')$value
tasmax.add_offset <- ncatt_get(nc.bcsd, varid='tasmax',
                               attname='add_offset')$value
tasmin.add_offset <- ncatt_get(nc.bcsd, varid='tasmin',
                               attname='add_offset')$value
wind.add_offset <- ncatt_get(nc.bcsd, varid='wind',
                             attname='add_offset')$value

################################################################################
# Calculate monthly climatologies of spatially aggregated observations

pr.aggregate.clim <- tasmax.aggregate.clim <- tasmin.aggregate.clim <-
    wind.aggregate.clim <- list()
for(mn in 1:12){
    pr.aggregate.clim.mn <- pr.aggregate.mn[(obs.time.mn[,1] %in%
                                resample.years) & (obs.time.mn[,2]==mn),]
    pr.aggregate.clim.mn <- colMeans(pr.aggregate.clim.mn)
    dim(pr.aggregate.clim.mn) <- dim(gcm.lons)
    pr.aggregate.clim[[mn]] <- pr.aggregate.clim.mn
    tasmax.aggregate.clim.mn <- tasmax.aggregate.mn[(obs.time.mn[,1] %in%
                                    resample.years) & (obs.time.mn[,2]==mn),]
    tasmax.aggregate.clim.mn <- colMeans(tasmax.aggregate.clim.mn)
    dim(tasmax.aggregate.clim.mn) <- dim(gcm.lons)
    tasmax.aggregate.clim[[mn]] <- tasmax.aggregate.clim.mn
    tasmin.aggregate.clim.mn <- tasmin.aggregate.mn[(obs.time.mn[,1] %in%
                                    resample.years) & (obs.time.mn[,2]==mn),]
    tasmin.aggregate.clim.mn <- colMeans(tasmin.aggregate.clim.mn)
    dim(tasmin.aggregate.clim.mn) <- dim(gcm.lons)
    tasmin.aggregate.clim[[mn]] <- tasmin.aggregate.clim.mn
    wind.aggregate.clim.mn <- wind.aggregate.mn[(obs.time.mn[,1] %in%
                                  resample.years) & (obs.time.mn[,2]==mn),]
    wind.aggregate.clim.mn <- colMeans(wind.aggregate.clim.mn)
    dim(wind.aggregate.clim.mn) <- dim(gcm.lons)
    wind.aggregate.clim[[mn]] <- wind.aggregate.clim.mn
}

################################################################################
# Dimensions for GCM downscaled outputs

for(end.day in c('31', '30', '29', '28')){
    gcm.time <- try(seq.Date(from=as.Date(paste(c(gcm.time.mn[1,1:2], '1'),
        collapse='-')), to=as.Date(paste(c(gcm.time.mn[nrow(gcm.time.mn),1:2],
        end.day), collapse='-')), by='days'), silent=TRUE)
    if(class(gcm.time) != 'try-error'){
        break
    }
}

ncvar_put(nc=nc.bcsd, varid='time', vals=seq_along(gcm.time)-1, start=1,
          count=length(gcm.time))
nc_sync(nc.bcsd)

n.lon <- nrow(obs.lons)
n.lat <- ncol(obs.lats)

################################################################################

set.seed(random.seed)
resample.time.mn <- obs.time.mn[obs.time.mn[,1] %in% resample.years,]
i.gcm <- 1
for(i in seq_along(gcm.time.mn[,1])){
    mn <- gcm.time.mn[i,2]
    cat(i.gcm, gcm.time.mn[i,1:2], ':')
    ##
    # Calculate and spatially disaggregate GCM scale/shift factors
    shift.pr.mn.z <- pr.gcm.mn[i,]/pr.aggregate.clim[[mn]]
    dim(shift.pr.mn.z) <- dim(gcm.lons)
    shift.pr.mn.grid <- list(x=gcm.lons[,1], y=gcm.lats[1,],
                             z=shift.pr.mn.z)
    shift.pr.mn <- interp.surface(shift.pr.mn.grid,
                                  cbind(c(obs.lons), c(obs.lats)))
    dim(shift.pr.mn) <- dim(obs.lons)

    shift.tasmax.mn.z <- tasmax.gcm.mn[i,] - tasmax.aggregate.clim[[mn]]
    dim(shift.tasmax.mn.z) <- dim(gcm.lons)
    shift.tasmax.mn.grid <- list(x=gcm.lons[,1], y=gcm.lats[1,],
                                 z=shift.tasmax.mn.z)
    shift.tasmax.mn <- interp.surface(shift.tasmax.mn.grid,
                                      cbind(c(obs.lons), c(obs.lats)))
    dim(shift.tasmax.mn) <- dim(obs.lons)

    shift.tasmin.mn.z <- tasmin.gcm.mn[i,] - tasmin.aggregate.clim[[mn]]
    dim(shift.tasmin.mn.z) <- dim(gcm.lons)
    shift.tasmin.mn.grid <- list(x=gcm.lons[,1], y=gcm.lats[1,],
                                 z=shift.tasmin.mn.z)
    shift.tasmin.mn <- interp.surface(shift.tasmin.mn.grid,
                                      cbind(c(obs.lons), c(obs.lats)))
    dim(shift.tasmin.mn) <- dim(obs.lons)

    shift.wind.mn.z <- wind.gcm.mn[i,]/wind.aggregate.clim[[mn]]
    dim(shift.wind.mn.z) <- dim(gcm.lons)
    shift.wind.mn.grid <- list(x=gcm.lons[,1], y=gcm.lats[1,],
                               z=shift.wind.mn.z)
    shift.wind.mn <- interp.surface(shift.wind.mn.grid,
                                    cbind(c(obs.lons), c(obs.lats)))
    dim(shift.wind.mn) <- dim(obs.lons)
    # Edge of domain correction due to bilinear interpolation from GCM back to
    # fine-scale grid
    if(i==1){
        edges <- which(is.na(shift.tasmin.mn) &
                       !is.na(tasmin.clim[[1]]))
        interior <- !is.na(shift.tasmin.mn)
        lonlat.interior <- cbind(obs.lons[interior], obs.lats[interior])
        nn.edges <- c()
        for(edge in edges){
            lonlat.edge <- c(obs.lons[edge], obs.lats[edge])
            nn <- which.min(rowSums(sweep(lonlat.interior, 2,
                            lonlat.edge, '-')^2))
            nn.edges <- c(nn.edges, which(interior)[nn])
        }
    }
    shift.pr.mn[edges] <- shift.pr.mn[nn.edges]
    shift.tasmax.mn[edges] <- shift.tasmax.mn[nn.edges]
    shift.tasmin.mn[edges] <- shift.tasmin.mn[nn.edges]
    shift.wind.mn[edges] <- shift.wind.mn[nn.edges]
    ##
    # Apply scale/shift factors to fine-scale climatology
    pr.disagg.mn <- shift.pr.mn*pr.clim[[mn]]
    if(shift.tavg){
        shift.tasmax.mn <- shift.tasmin.mn <-
            (shift.tasmax.mn + shift.tasmin.mn)/2
    }
    tasmax.disagg.mn <- shift.tasmax.mn + tasmax.clim[[mn]]
    tasmin.disagg.mn <- shift.tasmin.mn + tasmin.clim[[mn]]
    wind.disagg.mn <- shift.wind.mn*wind.clim[[mn]]
    ##
    # Temporally disaggregate by sampling month from observational record
    # Includes error checking on precipitation scaling factor and number of
    # wet days following Maurer and Hidalgo (2008)
    error.disagg <- TRUE
    pr.max.best <- Inf
    years.resample <- sample(which(resample.time.mn[,2]==mn))
    tries <- 1
    while(error.disagg){
        cat('*')
        i.resample <- years.resample[tries]
        yr.mn <- sapply(strsplit(names(i.resample), '-')[[1]], as.numeric)
        i.days <- which((obs.time[,1]==yr.mn[1]) & (obs.time[,2]==yr.mn[2]))
        # Correct for mismatch in number of days in month
        mn.char <- as.character(mn)
        if(nchar(mn.char)==1) mn.char <- paste('0', mn.char, sep='')
        n.mn.gcm <- sum(format(gcm.time, '%Y-%m')==
                        paste(gcm.time.mn[i,1], '-', mn.char, sep=''))
        i.days <- seq(i.days[1], length=n.mn.gcm)
        ##
        # Scale/shift daily precipitation values to match monthly means
        pr.obs <- ncvar_get(nc.obs, varid='pr', start=c(1, 1, i.days[1]),
                            count=c(n.lon, n.lat, length(i.days)))
        shift.pr.mn <- pr.disagg.mn/apply(pr.obs, c(1, 2), mean)
        pr.disagg <- sweep(pr.obs, c(1, 2), shift.pr.mn, '*')
        # Replace NaN values due to a dry month with zeros
        pr.disagg[which(is.nan(pr.disagg))] <- 0
        pr.disagg[which(pr.disagg < 0)] <- 0
        # Check for unrealistic precipitation scaling values
        pr.scale.resample <- apply(shift.pr.mn, c(1, 2), max)
        pr.nwet.resample <- apply(pr.disagg > 0, c(1, 2), sum)
        error.disagg <- sum((pr.scale.resample > pr.scale.max) &
                            (pr.nwet.resample <= pr.nwet.min), na.rm=TRUE) > 0
        error.disagg <- ((error.disagg) |
                         (max(pr.disagg, na.rm=TRUE)
                          > pr.prop.max*max(pr.max[[mn]], na.rm=TRUE)))
        if(!error.disagg){
            pr.max.best <- max(pr.disagg, na.rm=TRUE)
        } else{
            tries <- tries + 1
            if(max(pr.disagg, na.rm=TRUE) < pr.max.best){
                pr.max.best <- max(pr.disagg, na.rm=TRUE)
                pr.disagg.best <- pr.disagg
            }
        }
        if(tries > length(years.resample)){
            cat(' failed to find an analog month ')
            pr.disagg <- pr.disagg.best
            error.disagg <- FALSE
        }
    }
    # Scale/shift daily values to match monthly means
    tasmax.obs <- ncvar_get(nc.obs, varid='tasmax', start=c(1, 1, i.days[1]),
                            count=c(n.lon, n.lat, length(i.days)))
    tasmin.obs <- ncvar_get(nc.obs, varid='tasmin', start=c(1, 1, i.days[1]),
                            count=c(n.lon, n.lat, length(i.days)))
    shift.tasmax.mn <- tasmax.disagg.mn - apply(tasmax.obs, c(1, 2), mean)
    shift.tasmin.mn <- tasmin.disagg.mn - apply(tasmin.obs, c(1, 2), mean)
    if(shift.tavg){
        shift.tasmin.mn <- shift.tasmax.mn <-
            (shift.tasmin.mn + shift.tasmax.mn)/2
    }
    tasmax.disagg <- sweep(tasmax.obs, c(1, 2), shift.tasmax.mn, '+')
    tasmin.disagg <- sweep(tasmin.obs, c(1, 2), shift.tasmin.mn, '+')
    # Correct temperature reversals if handling tmin/tmax separately
    if(!shift.tavg){
        reversal <- which(tasmin.disagg > tasmax.disagg)
        tasmax.disagg.tmp <- tasmax.disagg
        tasmin.disagg.tmp <- tasmin.disagg
        tasmin.disagg[reversal] <- tasmax.disagg.tmp[reversal]
        tasmax.disagg[reversal] <- tasmin.disagg.tmp[reversal]
    }
    ##
    wind.obs <- ncvar_get(nc.obs, varid='wind', start=c(1, 1, i.days[1]),
                          count=c(n.lon, n.lat, length(i.days)))
    shift.wind.mn <- wind.disagg.mn/apply(wind.obs, c(1, 2), mean)
    wind.disagg <- sweep(wind.obs, c(1, 2), shift.wind.mn, '*')
    ##
    cat(':', yr.mn, max(pr.max[[mn]], na.rm=TRUE), pr.max.best, '\n')
    # Filter out overflow precipitation values
    pr.disagg[pr.disagg > pr.overflow] <- pr.overflow
     # Create packed data values
    pr.disagg <- round((pr.disagg - pr.add_offset)/pr.scale_factor)
    tasmax.disagg <- round((tasmax.disagg - tasmax.add_offset)/
                            tasmax.scale_factor)
    tasmin.disagg <- round((tasmin.disagg - tasmin.add_offset)/
                            tasmin.scale_factor)
    wind.disagg <- round((wind.disagg - wind.add_offset)/wind.scale_factor)
    # Missing values
    pr.disagg[is.na(pr.disagg)] <- pr.missing_value
    tasmax.disagg[is.na(tasmax.disagg)] <- tasmax.missing_value
    tasmin.disagg[is.na(tasmin.disagg)] <- tasmin.missing_value
    wind.disagg[is.na(wind.disagg)] <- wind.missing_value
    # Write packed data to NetCDF files
    ncvar_put(nc.bcsd, varid='pr', vals=pr.disagg,
              start=c(1, 1, i.gcm), count=c(n.lon, n.lat, length(i.days)))
    ncvar_put(nc.bcsd, varid='tasmax', vals=tasmax.disagg,
              start=c(1, 1, i.gcm), count=c(n.lon, n.lat, length(i.days)))
    ncvar_put(nc.bcsd, varid='tasmin', vals=tasmin.disagg,
              start=c(1, 1, i.gcm), count=c(n.lon, n.lat, length(i.days)))
    ncvar_put(nc.bcsd, varid='wind', vals=wind.disagg,
              start=c(1, 1, i.gcm), count=c(n.lon, n.lat, length(i.days)))
    i.gcm <- i.gcm + length(i.days)
    ##
    nc_sync(nc.bcsd)
}

nc_close(nc.obs)
nc_close(nc.bcsd)

##******************************************************************************
