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

compute.time.stats <- function(nc, start=NULL, end=NULL) {
  vals <- netcdf.calendar(nc, 'time')
  if (is.null(start)) {
    start <- vals[1]
  }
  if (is.null(end)) {
    end <- vals[length(vals)]
  }
  t0 <- as.PCICt(start, cal=attr(vals, 'cal'))
  tn <- as.PCICt(end, cal=attr(vals, 'cal'))
  i <- vals >= t0 & vals <= tn
  vals <- vals[i]
  list(vals=vals,
       i=i,
       t0=min(which(i)),
       tn=max(which(i)),
       n=length(vals)
       )
}
