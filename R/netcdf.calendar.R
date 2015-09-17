##******************************************************************************
# Convert NetCDF time dimension to a YMDH calendar
# Alex Cannon (acannon@uvic.ca)
##******************************************************************************

library(ncdf4)
library(PCICt)
library(RNetCDF)

netcdf.calendar <-
function(nc, time.variable='time', pcict=FALSE)
{
    time.calendar <- ncatt_get(nc, time.variable, 'calendar')$value
    if(time.calendar=='noleap') time.calendar <- '365_day'
    if(time.calendar==0) time.calendar <- 'gregorian'
    if(time.calendar=='standard') time.calendar <- 'gregorian'
    if(grepl('gregorian', time.calendar)){
        time.ymdh <- utcal.nc(ncatt_get(nc, time.variable, 'units')$value,
                              ncvar_get(nc, time.variable))[,1:4]
    } else{
        time.units <- ncatt_get(nc, time.variable, 'units')$value
        time.values <- ncvar_get(nc, time.variable)
        if(grepl('days', time.units)) time.values <- time.values*86400
        if(grepl('hours', time.units)) time.values <- time.values*3600
        origin.pcict <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                                 cal=time.calendar)
        time.values <- origin.pcict + time.values
        time.ymdh <- cbind(as.numeric(format(time.values, '%Y')),
                           as.numeric(format(time.values, '%m')),
                           as.numeric(format(time.values, '%d')),
                           as.numeric(format(time.values, '%H')))
        colnames(time.ymdh) <- c('year', 'month', 'day', 'hour')
    }
    if (pcict) {
        time.values
    } else {
        time.ymdh
    }
}
