##******************************************************************************
##******************************************************************************

##Modified to operate on Westgrid

library(ncdf4)
library(PCICt)
library(fields)
##******************************************************************************
################################################################################
##Calculate GCM anomalies

create.anoms <- function(var.name,infile,outfile,var.mon) {

  ##For mean subset 1951-2005
  file.copy(from=infile,to=outfile,overwrite=T)
  Sys.sleep(5)
  if (!file.exists(outfile))
    warning('Anomaly copy unsuccessful')
  
  nc <- nc_open(infile)
  var.data <- ncvar_get(nc,var.name)
  time.atts <- ncatt_get(nc,'time')
  time.calendar <- time.atts$calendar
  time.units <- time.atts$units
  time.values <- ncvar_get(nc,'time')
  origin.pcict <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                           cal=time.calendar)  
  var.dates <- origin.pcict + time.values*86400
  monthly.fac <- as.factor(format(var.dates,'%m'))
  
  ##Load full time series and take anomalies from this
  var.anoms <- var.data*0
  for(mn in 1:12) {
    print(mn)
    var.mean <- var.mon[mn,,]
    var.ix <- which(monthly.fac==sprintf('%02d',mn))
    mlen <- length(var.ix)
    for (i in 1:mlen) {
      ix <- var.ix[i]
      if (var.name=='pr')
        var.anoms[,,ix] <- var.data[,,ix]/var.mean
      if (grepl('tas',var.name))
        var.anoms[,,ix] <- var.data[,,ix] - var.mean
    }
  }
  nc_close(nc)
  ##Fix the edges for interpolations
  ncol <- dim(var.data)[2]
  for (j in 1:(ncol-1)) {
    ix <- is.na(var.anoms[,j,1])
    var.anoms[ix,j,] <- var.anoms[ix,(j+1),]    
  }
  
  if (var.name=='pr')
    var.anoms[is.na(var.anoms)] <- 0
  anc <- nc_open(anoms.file,write=TRUE)
  ncvar_put(anc,varid=var.name,vals=var.anoms)
  nc_close(anc)
}

##-------------------------------------------------------------

bccaq.anomalies <- function(var.name,gcm,gcm.file) {
  
  base.file <- gcm.file ##subset. file
  gnc <- nc_open(base.file)
  time.atts <- ncatt_get(gnc,'time')
  time.calendar <- time.atts$calendar
  time.units <- time.atts$units
  time.values <- ncvar_get(gnc,'time')
  origin.pcict <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                           cal=time.calendar)  
  var.dates <- origin.pcict + time.values*86400
  monthly.fac <- as.factor(format(var.dates,'%m'))
  monthly.ts.fac <- as.factor(format(var.dates,'%Y-%m'))
  mon.facs <- as.factor(format(as.Date(paste(levels(monthly.ts.fac),'-01',sep='')),'%m'))
  var.data <- ncvar_get(gnc,var.name)
  var.mon <- apply(var.data,c(1,2),function(x,fac){tapply(x,fac,mean,na.rm=T)},monthly.fac)
  if (var.name=='pr') {
    var.data[var.data <=0] <- NA    
    var.test <- apply(var.data,c(1,2),function(x,fac){tapply(x,fac,sum,na.rm=T)},monthly.ts.fac)
    var.mon <-  apply(var.test,c(2,3),function(x,fac){tapply(x,fac,mean,na.rm=T)},mon.facs)
  }
  nc_close(gnc)
  ##------------------------------------  
  create.anoms(var.name,infilebase.file,var.mon) ##var.mon is the monthly gcm climatologies

}

interp.bccaq <- function(var.name,gcm,interval,base.dir) {
  print(paste('Interpolate Anomalies: ',gcm,', ',var.name,', ',interval,sep=''))
  ##grid.file <- '/home/data/projects/rci/data/stat.downscaling/BCCAQ/bccaq_gcm_van_whistler_subset/prism.grid.txt'
  grid.file <- paste(base.dir,'prism.grid.txt',sep='')
  gcm.dir <- paste(base.dir,gcm,'/',sep='')  
  ##gcm.dir <- paste('/home/data/projects/rci/data/stat.downscaling/BCCAQ/bccaq_gcm_van_whistler_subset/',gcm,'/',sep='')
  anoms.files <- list.files(path=gcm.dir,pattern=paste(var.name,'_anoms_BCCAQ_',sep=''),full.name=TRUE)
  anoms.file <- anoms.files[grep(interval,anoms.files)]
  interp.file <- gsub(pattern='_anoms_BCCAQ_',replacement='_anoms_interp_BCCAQ_',anoms.file)
  system(paste('cdo -s remapbil,',grid.file,' ',anoms.file,' ',interp.file,sep=''))

  if (1==0) {
    ##Test of interp.grid
    tasmax.nc <- nc_open(anoms.file)
    lat <- ncvar_get(tasmax.nc, 'lat')
    lon <- ncvar_get(tasmax.nc, 'lon')
    lonlat <- list(lon=lon, lat=lat)
    tasmax.test <- ncvar_get(tasmax.nc,'tasmax',start=c(1,1,1),count=c(-1,-1,1))
    tx.grid <- list(x=lonlat$lon, y=lonlat$lat, z=tasmax.test)
    
    obs.lat <- seq(from=48.91667,by=0.0083333,length.out=180)
    obs.lon <- seq(from=-123.5083,by=0.0083333,length.out=339)
    obs.lats <- c(matrix(obs.lat, ncol=length(obs.lat),
                         nrow=length(obs.lon), byrow=TRUE))
    obs.lons <- c(matrix(obs.lon, ncol=length(obs.lat),
                         nrow=length(obs.lon)))
    obs.coords <- cbind(obs.lons, obs.lats)
    
    bilin.tx <- interp.surface(tx.grid, obs.coords)
    browser()
  }
  
}

daily.prism.scale <- function(var.name,gcm,interval,base.dir) {
  print(paste('Daily PRISM: ',gcm,', ',var.name,', ',interval,sep=''))
  ##base.dir <- '/home/data/projects/rci/data/stat.downscaling/BCCAQ/bccaq_gcm_van_whistler_subset/'
  gcm.dir <- paste(base.dir,gcm,'/',sep='')
  bccaq.files <- list.files(path=gcm.dir,pattern=paste(var.name,'_anoms_interp_BCCAQ_',sep=''),full.name=TRUE)
  bccaq.file <- bccaq.files[grep(interval,bccaq.files)]
  prism.var <- switch(var.name,
                      pr='pr',
                      tasmax='tmax',
                      tasmin='tmin')
  prism.file <- paste(base.dir,'PRISM/',prism.var,'_monClim_PRISM_historical_van_whistler_1971-2000.nc',sep='')
  ##prism.file <- paste('/home/data/projects/rci/data/stat.downscaling/BCCAQ/bccaq_gcm_van_whistler_subset/PRISM/',
  ##                    prism.var,'_monClim_PRISM_historical_van_whistler_1971-2000.nc',sep='')

  adjusted.file <- gsub(pattern='_anoms_interp_BCCAQ_',replacement='_gcm_prism_BCCAQ_',bccaq.file)
  file.copy(from=bccaq.file,to=adjusted.file,overwrite=T)  
  ##PRISM climatologies
  pnc <- nc_open(prism.file)
  prism.clim <- ncvar_get(pnc,prism.var)
  nc_close(pnc)
  
  bnc <- nc_open(bccaq.file)
  anc <- nc_open(adjusted.file,write=TRUE)
  
  time.atts <- ncatt_get(bnc,'time')
  time.calendar <- time.atts$calendar
  time.units <- time.atts$units
  time.values <- ncvar_get(bnc,'time')
  origin.pcict <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                           cal=time.calendar)  
  var.dates <- origin.pcict + time.values*86400
  monthly.fac <- as.factor(format(var.dates,'%m'))

##Flags for the SW corner
  flag.file <- list.files(path=paste(base.dir,gcm,sep=''),pattern='tasmax_gcm_prism_BCCAQ_',full.name=TRUE)[1]
  fnc <- nc_open(flag.file,write=TRUE)
  flag.data <- ncvar_get(fnc,'tasmax',start=c(1,1,1),count=c(-1,-1,1))
  flags <- is.na(flag.data)
  nc_close(fnc)
  
  ##Break up the data into pieces to be manageable
  tlen <- length(time.values)
  sqnce <- seq(0,100000,by=4000)
  sx <- findInterval(tlen,sqnce)
  itvls <- c(sqnce[1:sx],tlen-sqnce[sx]+sqnce[sx])

  nlon <- bnc$dim$lon$len
  nlat <- bnc$dim$lat$len
  addon <- matrix(NA,nrow=339,ncol=9)
  
  for (s in 1:sx) {
    st <- itvls[s] + 1
    en <- itvls[s+1]
    len <- en-st+1
    print(paste(st,' to ',en,' of ',tlen,sep=''))
    var.data <- ncvar_get(bnc,var.name,start=c(1,1,st),count=c(-1,-1,len))

    var.adjust <- var.data*0  
    fac.sub <- monthly.fac[st:en]
    for(mn in 1:12) {
      print(mn)
      prism.mean <- cbind(addon,prism.clim[,,mn])
      
      var.ix <- which(fac.sub==sprintf('%02d',mn))
      mlen <- length(var.ix)
      for (i in 1:mlen) {
        ix <- var.ix[i]
        var.sub <- var.adjust[,,ix]
        if (var.name=='pr')
          var.sub <- var.data[,,ix]*prism.mean
          ##var.adjust[,10:nlat,ix] <- var.data[,10:nlat,ix]*prism.mean
        if (grepl('tas',var.name))
          var.sub <- var.data[,,ix] + prism.mean
          ##var.adjust[,10:nlat,ix] <- var.data[,10:nlat,ix] + prism.mean
        var.sub[flags] <- NA        
        var.adjust[,,ix] <- var.sub
      }##Loop over indices    
    }##Loop over Months
    ##Flag the NA values
    ##browser()
    ##if (var.name=='pr') {
    ##  var.adjust[flags] <- NA
    ##}
    ##var.adjust[,1:9,] <- NA

    ncvar_put(anc,varid=var.name,vals=var.adjust,start=c(1,1,st),count=c(nlon,nlat,len))
  }##Loop over data pieces
  rm(var.adjust)
  nc_close(bnc)
  nc_close(anc)
}
################################################################################
##******************************************************************************

##First 3 stages of BCCI to get climate imprint

code.dir <- '/home/ssobie/stat.downscaling/code/QPQM/BCCI/'

config <- paste(code.dir,'BCCI.set.config',sep='')
print(readLines(config))
source(config)

var.name <- 'pr'
var.title <- 'precip'
bccaq.anomalies(var.name,gcm,base.dir)
interp.bccaq(var.name,gcm,'1951-2000',base.dir)
daily.prism.scale(var.name,gcm,'1951-2000',base.dir)




