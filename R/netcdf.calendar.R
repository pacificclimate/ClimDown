# Convert NetCDF time dimension to a PCICt vector

netcdf.calendar <-
function(nc, time.variable='time') {
    time.calendar <- ncatt_get(nc, time.variable, 'calendar')$value
    time.units <- ncatt_get(nc, time.variable, 'units')$value
    time.values <- ncvar_get(nc, time.variable)

    origin.pcict <- as.PCICt(strsplit(time.units, ' ')[[1]][3],
                             cal=time.calendar)

    if (time.calendar == 'noleap') {
        time.calendar <- '365_day'
    }
    else if (time.calendar == 0 || time.calendar == 'standard') {
        time.calendar <- 'gregorian'
    }
    if (grepl('days', time.units)) {
        time.values <- time.values * 86400
    } else if(grepl('hours', time.units)) {
        time.values <- time.values * 3600
    }
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
  i <- compute.time.overlap(vals, t0, tn)
  vals <- vals[i]
  list(vals=vals,
       i=i,
       t0=min(which(i)),
       tn=max(which(i)),
       n=length(vals)
       )
}

compute.time.overlap <- function(timevals, t0, tn, error=TRUE) {
    ti <- timevals >= t0 & timevals <= tn
    if (! any(ti)) {
        d <- format(c(t0, tn, range(timevals)), '%Y-%m-%d')
        msg <- paste(
            c("The configured calibration period (", d[1:2],
              ") does not overlap with the GCM time period (", d[3:4], ")"),
            collapse=' '
            )
        if (error) {
            stop(msg)
        } else {
            warning(msg)
        }
    }
    ti
}
