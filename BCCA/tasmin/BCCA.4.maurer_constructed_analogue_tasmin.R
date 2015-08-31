##******************************************************************************
# Bias Corrected Constructed Analogue (BCCA) downscaling algorithm
# Alex Cannon (acannon@uvic.ca)
# *All wind related lines of code - commented out by Arelia Werner (wernera@uvic.ca)
##******************************************************************************
# Constructed analogues for each bias-corrected GCM day
##******************************************************************************
ptm <- proc.time()
print('Starting')

code.dir <- '/home/ssobie/stat.downscaling/code/QPQM/BCCA/'

source(paste(code.dir,'DQM.R',sep=''))
library(ncdf4)

args <- commandArgs(trailingOnly=TRUE)
for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
}

config <- paste(code.dir,'BCCA.set.config',sep='')
print(readLines(config))
source(config)

nc <- nc_open(nc.obs.file)

system(paste('cp ', template.file, '_tasmin_only.nc ', output.dir, output.file,
       output.suffix, '_tasmin.nc', sep=''))

nc.bcca <- nc_open(paste(output.dir, output.file, output.suffix, '_tasmin.nc',
                   sep=''), write=TRUE)

ncatt_put(nc.bcca, varid=0, attname='title', attval=output.nc.title)
ncatt_put(nc.bcca, varid=0, attname='institution',
          attval=output.nc.institution)
ncatt_put(nc.bcca, varid=0, attname='source', attval=output.nc.source)
ncatt_put(nc.bcca, varid=0, attname='input_data', attval=output.nc.input_data)
ncatt_put(nc.bcca, varid=0, attname='reference', attval=output.nc.reference)
ncatt_put(nc.bcca, varid=0, attname='project_id', attval=output.nc.project_id)
ncatt_put(nc.bcca, varid=0, attname='experiment_id',
          attval=output.nc.experiment_id)
ncatt_put(nc.bcca, varid=0, attname='version', attval=output.nc.version)
ncatt_put(nc.bcca, varid=0, attname='version_comment',
          attval=output.nc.version_comment)
ncatt_put(nc.bcca, varid=0, attname='contact1', attval=output.nc.contact1)
ncatt_put(nc.bcca, varid=0, attname='contact2', attval=output.nc.contact2)
ncatt_put(nc.bcca, varid=0, attname='contact3', attval=output.nc.contact3)
ncatt_put(nc.bcca, varid=0, attname='history', attval=output.nc.history)
ncatt_put(nc.bcca, varid='time', attname='units',
          attval=output.nc.time.units)
ncatt_put(nc.bcca, varid='time', attname='calendar',
          attval=output.nc.calendar)
nc_sync(nc.bcca)

load(paste(output.dir, 'gcm.lons', output.suffix, '.RData', sep=''))
load(paste(output.dir, 'gcm.lats', output.suffix, '.RData', sep=''))
load(paste(output.dir, 'obs.lons', output.suffix, '.RData', sep=''))
load(paste(output.dir, 'obs.lats', output.suffix, '.RData', sep=''))
load(paste(output.dir, 'obs.time', output.suffix, '.RData', sep=''))
load(paste(output.dir, 'tasmin.aggregate', output.suffix, '.RData', sep=''))

load(paste(output.dir, 'tasmin.gcm.time', output.suffix, '.RData', sep=''))
load(paste(output.dir, 'tasmin.raw.time', output.suffix, '.RData', sep=''))
load(paste(output.dir, 'tasmin.gcm_bc', output.suffix, '.RData', sep=''))

##******************************************************************************

tasmin.missing_value <- ncatt_get(nc.bcca, varid='tasmin',
                              attname='missing_value')$value

tasmin.scale_factor <- ncatt_get(nc.bcca, varid='tasmin',
                             attname='scale_factor')$value

tasmin.add_offset <- ncatt_get(nc.bcca, varid='tasmin',
                           attname='add_offset')$value

################################################################################

n.lon <- nrow(obs.lons)
n.lat <- ncol(obs.lats)

obs.Date <- as.Date(paste(1996, obs.time[,2], obs.time[,3], sep='-'))
gcm.D1 <- as.Date(paste(1995, tasmin.gcm.time[,2], tasmin.gcm.time[,3], sep='-'))
gcm.D2 <- as.Date(paste(1996, tasmin.gcm.time[,2], tasmin.gcm.time[,3], sep='-'))
gcm.D3 <- as.Date(paste(1997, tasmin.gcm.time[,2], tasmin.gcm.time[,3], sep='-'))

if (gcm == 'hadcm3') {
  gcm.D1[is.na(gcm.D1)] <- as.Date('1995-02-28')
  gcm.D2[is.na(gcm.D2)] <- as.Date('1996-02-28')
  gcm.D3[is.na(gcm.D3)] <- as.Date('1997-02-28')
}


na.mask <- !is.na(tasmin.gcm[1,])
for(i in seq_along(tasmin.gcm.time[,1])){
    ncvar_put(nc=nc.bcca, varid='time', vals=tasmin.raw.time[i], start=i, count=1)
    # Develop library of observed days within +/- delta.days of the
    # GCM simulated day
    alib <- which(((obs.Date <= (gcm.D1[i] + delta.days)) &
                   (obs.Date >= (gcm.D1[i] - delta.days))) |
                  ((obs.Date <= (gcm.D2[i] + delta.days)) &
                   (obs.Date >= (gcm.D2[i] - delta.days))) |
                  ((obs.Date <= (gcm.D3[i] + delta.days)) &
                   (obs.Date >= (gcm.D3[i] - delta.days))))
    alib <- alib[obs.time[alib,1] %in% obs.ca.years]
    ## Maximum Temperature
    # Find the n.analogue closest observations from the library
    tasmin.gcm.i <- tasmin.gcm[i,na.mask]
    tasmin.agg.alib <- tasmin.aggregate[alib,na.mask]
    analogues <- which(rank(rowSums(sweep(tasmin.agg.alib, 2, tasmin.gcm.i, '-')^2),
                       ties.method='random') %in% 1:n.analogue)
    # Constructed analogue weights
    tasmin.agg.alib <- jitter(tasmin.agg.alib[analogues,])
    # tasmin.agg.alib <- tasmin.agg.alib[analogues,]
    tasmin.Q <- tasmin.agg.alib %*% t(tasmin.agg.alib)
    tasmin.ridge <- tol*mean(diag(tasmin.Q))
    tasmin.ridge <- diag(n.analogue)*tasmin.ridge
    tasmin.weights <- (solve(tasmin.Q + tasmin.ridge) %*%
                   tasmin.agg.alib) %*% as.matrix(tasmin.gcm.i)
    tasmin.analogue <- 0
    for(j in 1:n.analogue){
        tasmin.analogue.j <- ncvar_get(nc=nc, varid='tasmin',
                                   start=c(1, 1, alib[analogues[j]]),
                                   count=c(n.lon, n.lat, 1))
        tasmin.analogue <- tasmin.analogue + tasmin.weights[j]*(tasmin.analogue.j)
    }
    cat('*')

    ##
    # Create packed data values
    tasmin.analogue <- round((tasmin.analogue - tasmin.add_offset)/tasmin.scale_factor)
    ##
    # Missing values
    tasmin.analogue[is.na(tasmin.analogue)] <- tasmin.missing_value

    # Write packed data to NetCDF files
    ncvar_put(nc.bcca, varid='tasmin', vals=tasmin.analogue,
              start=c(1, 1, i), count=c(n.lon, n.lat, 1))
    nc_sync(nc.bcca)
    ##
}

################################################################################

nc_close(nc)
nc_close(nc.bcca)

print('Elapsed Time')
print(proc.time() - ptm)