##******************************************************************************
# Convert NetCDF time dimension to a YMDH calendar
# Alex Cannon (acannon@uvic.ca)
##******************************************************************************

#library(RNetCDF)

netcdf.calendar <-
function(nc, time.variable='time')
{
    time.calendar <- ncatt_get(nc, time.variable, 'calendar')$value
    time.units <- ncatt_get(nc, time.variable, 'units')$value
    time.values <- ncvar_get(nc, time.variable)

    origin.pcict <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                             cal=time.calendar)

    if(time.calendar=='noleap') time.calendar <- '365_day'
    if(time.calendar==0) time.calendar <- 'gregorian'
    if(time.calendar=='standard') time.calendar <- 'gregorian'
    if(grepl('days', time.units)) time.values <- time.values*86400
    if(grepl('hours', time.units)) time.values <- time.values*3600
    origin.pcict + time.values
}
