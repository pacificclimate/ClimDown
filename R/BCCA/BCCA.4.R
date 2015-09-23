##******************************************************************************
# Bias Corrected Constructed Analogue (BCCA) downscaling algorithm
# Alex Cannon (acannon@uvic.ca)
# *All wind related lines of code - commented out by Arelia Werner (wernera@uvic.ca)
##******************************************************************************
# Constructed analogues for each bias-corrected GCM day
##******************************************************************************

library(ncdf4)
code.dir <- Sys.getenv('CODE_DIR')
source(paste(code.dir, 'netcdf.calendar.R', sep='/'))
source(paste(code.dir, 'BCCA/BCCA.R', sep='/'))

# "--args gcm.file='${gcm.file}' obs.file='${obs.file}' output.file='${output.file}' varid='${varid}'"
args <- commandArgs(trailingOnly=TRUE)
for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
}

nc.obs <- nc_open(obs.file)
nc.gcm <- nc_open(gcm.file)

#system(paste('cp ', template.file, '_tasmax_only.nc ', output.dir, output.file,
#       output.suffix, '_tasmax.nc', sep=''))

nc.bcca <- nc_open(output.file, write=TRUE)

## ncatt_put(nc.bcca, varid=0, attname='title', attval=output.nc.title)
## ncatt_put(nc.bcca, varid=0, attname='institution',
##           attval=output.nc.institution)
## ncatt_put(nc.bcca, varid=0, attname='source', attval=output.nc.source)
## ncatt_put(nc.bcca, varid=0, attname='input_data', attval=output.nc.input_data)
## ncatt_put(nc.bcca, varid=0, attname='reference', attval=output.nc.reference)
## ncatt_put(nc.bcca, varid=0, attname='project_id', attval=output.nc.project_id)
## ncatt_put(nc.bcca, varid=0, attname='experiment_id',
##           attval=output.nc.experiment_id)
## ncatt_put(nc.bcca, varid=0, attname='version', attval=output.nc.version)
## ncatt_put(nc.bcca, varid=0, attname='version_comment',
##           attval=output.nc.version_comment)
## ncatt_put(nc.bcca, varid=0, attname='contact1', attval=output.nc.contact1)
## ncatt_put(nc.bcca, varid=0, attname='contact2', attval=output.nc.contact2)
## ncatt_put(nc.bcca, varid=0, attname='contact3', attval=output.nc.contact3)
## ncatt_put(nc.bcca, varid=0, attname='history', attval=output.nc.history)
## ncatt_put(nc.bcca, varid='time', attname='units',
##           attval=output.nc.time.units)
## ncatt_put(nc.bcca, varid='time', attname='calendar',
##           attval=output.nc.calendar)
## nc_sync(nc.bcca)

## load(paste(output.dir, 'gcm.lons', output.suffix, '.RData', sep=''))
## load(paste(output.dir, 'gcm.lats', output.suffix, '.RData', sep=''))
## load(paste(output.dir, 'obs.lons', output.suffix, '.RData', sep=''))
## load(paste(output.dir, 'obs.lats', output.suffix, '.RData', sep=''))
## load(paste(output.dir, 'obs.time', output.suffix, '.RData', sep=''))
## load(paste(output.dir, 'tasmax.aggregate', output.suffix, '.RData', sep=''))

## load(paste(output.dir, 'tasmax.gcm.time', output.suffix, '.RData', sep=''))
## load(paste(output.dir, 'tasmax.raw.time', output.suffix, '.RData', sep=''))
## load(paste(output.dir, 'tasmax.gcm_bc', output.suffix, '.RData', sep=''))

##******************************************************************************

## tasmax.missing_value <- ncatt_get(nc.bcca, varid='tasmax',
##                               attname='missing_value')$value

## tasmax.scale_factor <- ncatt_get(nc.bcca, varid='tasmax',
##                              attname='scale_factor')$value

## tasmax.add_offset <- ncatt_get(nc.bcca, varid='tasmax',
##                            attname='add_offset')$value

################################################################################

obs.time <- netcdf.calendar(nc.obs)
gcm.time <- netcdf.calendar(nc.gcm)

obs.Date <- as.Date(paste(1996, obs.time[,2], obs.time[,3], sep='-'))
gcm.D1 <- as.Date(paste(1995, gcm.time[,2], gcm.time[,3], sep='-'))
gcm.D2 <- as.Date(paste(1996, gcm.time[,2], gcm.time[,3], sep='-'))
gcm.D3 <- as.Date(paste(1997, gcm.time[,2], gcm.time[,3], sep='-'))

## if (gcm == 'hadcm3') {
##   gcm.D1[is.na(gcm.D1)] <- as.Date('1995-02-28')
##   gcm.D2[is.na(gcm.D2)] <- as.Date('1996-02-28')
##   gcm.D3[is.na(gcm.D3)] <- as.Date('1997-02-28')
## }

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


# returns 30 indices of the timestep for the closest analog and their corresponding weights
find.analogues <- function(gcm, agg) {
    na.mask <- !is.na(gcm)
    alib <- which(((obs.Date <= (gcm.D1[i] + delta.days)) &
                       (obs.Date >= (gcm.D1[i] - delta.days))) |
                           ((obs.Date <= (gcm.D2[i] + delta.days)) &
                                (obs.Date >= (gcm.D2[i] - delta.days))) |
                                    ((obs.Date <= (gcm.D3[i] + delta.days)) &
                                         (obs.Date >= (gcm.D3[i] - delta.days))))
    alib <- alib[obs.time[alib,1] %in% obs.ca.years] ## FIXME: alib is only be defined for each julian day in any year... it is 50x redundant as defined

    # Find the n.analogue closest observations from the library
    gcm.i <- gcm[i, na.mask] # A single time step of GCM values
    # (obs years * (delta days * 2 + 1)) x cells
    tasmax.agg.alib <- tasmax.aggregate[alib, na.mask] # FIXME: this is the same for each julian day of any year
    # substract the GCM at this time step from the aggregated obs *for every library time value*
    # square that difference and then find the 30 lowest differences
    # returns n analogues for this particular GCM timestep
    analogues <- which(rank(rowSums(sweep(tasmax.agg.alib, 2, tasmax.gcm.i, '-')^2), # FIXME: colsums is faster
                       ties.method='random') %in% 1:n.analogue) # FIXME, replace which(rank) with sort()[1:30]
    # Constructed analogue weights
    weights <- construct.analogue.weights(tasmax.agg.alib[analogues,], tasmax.gcm.i)
    list(analogues=analogues, weights=weigths)
}

n.analogues <- 30
delta.days <- 45
obs.ca.years <- 1951:2005
tol <- 0.1 ##0.001
expon <- 0.5

gcm <- ncvar_get(nc.gcm, varid)-273.15 # FIXME: This is tas, but make this call udunits
gcm.time <- netcdf.calendar(nc.gcm, 'time', pcict=TRUE)
#gcm <- sweep(gcm, 2, na.mask, '*') # Pretty sure that this is unnecessary
gcm <- round(gcm, 3) # FIXME: This depends on the variable

# Read the aggregated obs
print(paste('Starting spatial aggregation', obs.file, gcm.file, varid))
ptm <- proc.time()
aggd.obs <- create.aggregates(obs.file, gcm.file, varid)
print('Elapsed time')
print(proc.time() - ptm)

print(paste('Starting bias correction', obs.file, gcm.file, varid))
ptm <- proc.time()
obs.time <- netcdf.calendar(nc.obs, 'time', pcict=TRUE)
gcm <- bias.correct.dqm(gcm, aggd.obs, obs.time, gcm.time, detrend=FALSE)
print('Elapsed time')
print(proc.time()-ptm)


print(paste('Starting analaogue construction', obs.file, gcm.file, varid))
ptm <- proc.time()
na.mask <- !is.na(aggd.obs[,,1])
Rprof()

obs.time <- netcdf.calendar(nc.obs)
for(i in seq_along(gcm.time)) {
    ncvar_put(nc=nc.bcca, varid='time', vals=i, start=i, count=1) # FIXME: do this all in one write outside of the loop
    # Develop library of observed days within +/- delta.days of the
    # GCM simulated day
    alib <- which(((obs.Date <= (gcm.D1[i] + delta.days)) &
                   (obs.Date >= (gcm.D1[i] - delta.days))) |
                  ((obs.Date <= (gcm.D2[i] + delta.days)) &
                   (obs.Date >= (gcm.D2[i] - delta.days))) |
                  ((obs.Date <= (gcm.D3[i] + delta.days)) &
                   (obs.Date >= (gcm.D3[i] - delta.days))))
    alib <- alib[obs.time[alib,1] %in% obs.ca.years] ## FIXME: alib is only be defined for each julian day in any year... it is 50x redundant as defined
    ## Maximum Temperature
    # Find the n.analogue closest observations from the library
    gcm.i <- gcm[,,i] # A single time step of GCM values
    # (obs years * (delta days * 2 + 1)) x cells
    agg.alib <- aggd.obs[,,alib] # FIXME: this is the same for each julian day of any year

    # substract the GCM at this time step from the aggregated obs *for every library time value*
    # square that difference and then find the 30 lowest differences
    # returns n analogues for this particular GCM timestep
    diffs <- sweep(agg.alib, 1:2, gcm.i, '-')^2
    # FIXME, replace which(rank) with sort()[1:30]
    analogues <- which(rank(apply(diffs, 3, sum, na.rm=T), ties.method='random') %in% 1:n.analogues)

    # Constructed analogue weights
    obs.at.analogues <- t(matrix(agg.alib[,,analogues][na.mask], ncol=n.analogues))
    weights <- construct.analogue.weights(obs.at.analogues, gcm.i[na.mask])

    # FIXME: This can easily be sum(mapply(weights, analog indices))
    analogue <- 0
    for(j in 1:n.analogues){
        analogue.j <- ncvar_get(nc=nc.obs, varid='tasmax',
                                start=c(1, 1, alib[analogues[j]]),
                                count=c(-1, -1, 1))
        analogue <- analogue + weights[j]*(analogue.j)
    }
    cat('*')
    ##
    # Create packed data values
    tasmax.add_offset <- 0
    tasmax.scale_factor <- 1.0
    tasmax.missing_value <- NA
    analogue <- round((analogue - tasmax.add_offset)/tasmax.scale_factor) # FIXME: Don't bother packing
    ##
    # Missing values
    analogue[is.na(analogue)] <- tasmax.missing_value # FIXME: Don't bother. The library should do this

    # Write packed data to NetCDF files
    ncvar_put(nc.bcca, varid=varid, vals=analogue,
              start=c(1, 1, i), count=c(-1, -1, 1))
    nc_sync(nc.bcca)
    ##
}
print('Elapsed time')
print(proc.time()-ptm)

Rprof(NULL)
summaryRprof()

################################################################################

nc_close(nc.gcm)
nc_close(nc.obs)
nc_close(nc.bcca)
