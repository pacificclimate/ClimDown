##******************************************************************************
# Bias Corrected Climate Imprint (BCCI) downscaling algorithm
# Alex Cannon (acannon@uvic.ca)
# *All wind related lines of code - commented out by Arelia Werner (wernera@uvic.ca)
##******************************************************************************
# Read GCM grid and calculate daily anomalies with respect to the GCM monthly
# climatologies
##******************************************************************************
##Modified to only analyze precipitation data
##******************************************************************************

ptm <- proc.time()

library(ncdf4)

code.dir <- '/home/ssobie/stat.downscaling/code/QPQM/BCCI/'

source(paste(code.dir,'netcdf.calendar.R',sep=''))

mail <- function(address, subject, message) {
  system(paste("echo '", message,
               "' | mail -s '", subject,
               "' ", address, sep=""))
}

config <- commandArgs(trailingOnly=TRUE)
if(length(config)==0) config <- paste(code.dir,'BCCI.set.config',sep='')
print(readLines(config))

  source(config)
  
  pr.nc <- nc_open(pr.nc.file)
  
  ##******************************************************************************
  ## Read GCM data
  
  if (!file.exists(output.dir))
     dir.create(output.dir,recursive=T)

  pr <- ncvar_get(pr.nc, pr.var)*86400
  pr <- round(pr, 3)

  lat <- ncvar_get(pr.nc, 'lat')
  lon <- ncvar_get(pr.nc, 'lon')
  lonlat <- list(lon=lon, lat=lat)
  save(lonlat, file=paste(output.dir, 'lonlat.gcm', output.suffix, '.RData',
                 sep=''))

  raw.time <- ncvar_get(pr.nc,'time')
  pr.dates <- netcdf.calendar(pr.nc)
  pr.clim.indices <- pr.dates[,1] %in% clim.years

  save(raw.time, file=paste(output.dir, 'pr.raw.time', output.suffix, '.RData',
                sep=''))

  dates <- pr.dates
  save(dates, file=paste(output.dir, 'pr.dates.gcm', output.suffix, '.RData',
                sep=''))

  nc_close(pr.nc)
  rm(raw.time)
  
################################################################################
# Calculate monthly GCM climatologies

  pr.clim <- list()

  for(mn in 1:12){
    pr.mn.clim.indices <- pr.clim.indices & (pr.dates[,2]==mn)
    pr.clim[[mn]] <- apply(pr[,,pr.mn.clim.indices], c(1,2), mean)
  }

################################################################################
# Calculate daily anomalies with respect to monthly GCM climatologies

  for(i in seq_along(pr.dates[,1])){
    cat(pr.dates[i,], '\n')
    mn <- pr.dates[i,2]
    pr[,,i] <- pr[,,i]/pr.clim[[mn]]
  }
  save(pr, file=paste(output.dir, 'pr.gcm', output.suffix, '.RData',
             sep=''))
  rm(pr)
  rm(pr.dates)
  rm(pr.clim)

print('Elapsed time')
print(proc.time() - ptm)

################################################################################
  ##Using tasmax dates for now
##mail("ssobie@uvic.ca","BCCI.2 finished",paste("The BCCI.2 script has finished with ",gcm,sep=""))


rm(list=ls())
##******************************************************************************
##}
