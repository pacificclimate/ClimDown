#!/usr/bin/Rscript

rm(list=ls())
library(ncdf4)

##******************************************************************************

forcings.file <- 'forcings_new_symap_BC_22AUG2011_1950-2006.nc'
template.file <- 'template.nc'

time.units <- 'days since 1948-01-01 00:00:00'
missing_value <- -2^(16-1)

pr.scale_factor <- 0.025
tasmax.scale_factor <- 0.01
tasmin.scale_factor <- 0.01
wind.scale_factor <- 0.01

pr.add_offset <- 750
tasmax.add_offset <- 0
tasmin.add_offset <- 0
wind.add_offset <- 0

nc.title <- 'Template file for spatial downscaling to VIC grid'
nc.institution <- 'Pacific Climate Impacts Consortium (PCIC), Victoria, BC, www.pacificclimate.org'
nc.source <- ''
nc.input_data <- ''
nc.reference <- ''
nc.project_id <- ''
nc.experiment_id <- ''
nc.version <- ''
nc.version_comment <- ''
nc.contact1 <- 'Alex Cannon'
nc.contact2 <- 'acannon@uvic.ca'
nc.contact3 <- ''
nc.history <- ''

##******************************************************************************

system(paste('rm -rf', template.file))
system(paste('ncks -H -d time,0', forcings.file, template.file))
nc <- nc_open(template.file, write=TRUE)

ncatt_put(nc, varid='time', attname='units', attval=time.units)
for(varid in c('pr', 'tasmax', 'tasmin', 'wind'))
    ncatt_put(nc, varid=varid, attname='missing_value', attval=missing_value)

ncatt_put(nc, varid='pr', attname='scale_factor',
          attval=pr.scale_factor, prec='float')
ncatt_put(nc, varid='tasmax', attname='scale_factor',
          attval=tasmax.scale_factor, prec='float')
ncatt_put(nc, varid='tasmin', attname='scale_factor',
          attval=tasmax.scale_factor, prec='float')
ncatt_put(nc, varid='wind', attname='scale_factor',
          attval=wind.scale_factor, prec='float')

ncatt_put(nc, varid='pr', attname='add_offset',
          attval=pr.add_offset, prec='float')
ncatt_put(nc, varid='tasmax', attname='add_offset',
          attval=tasmax.add_offset, prec='float')
ncatt_put(nc, varid='tasmin', attname='add_offset',
          attval=tasmax.add_offset, prec='float')
ncatt_put(nc, varid='wind', attname='add_offset',
          attval=wind.add_offset, prec='float')

ncatt_put(nc, varid=0, attname='title', attval=nc.title)
ncatt_put(nc, varid=0, attname='institution', attval=nc.institution)
ncatt_put(nc, varid=0, attname='source', attval=nc.source)
ncatt_put(nc, varid=0, attname='input_data', attval=nc.input_data)
ncatt_put(nc, varid=0, attname='reference', attval=nc.reference)
ncatt_put(nc, varid=0, attname='project_id', attval=nc.project_id)
ncatt_put(nc, varid=0, attname='experiment_id', attval=nc.experiment_id)
ncatt_put(nc, varid=0, attname='version', attval=nc.version)
ncatt_put(nc, varid=0, attname='version_comment', attval=nc.version_comment)
ncatt_put(nc, varid=0, attname='contact1', attval=nc.contact1)
ncatt_put(nc, varid=0, attname='contact2', attval=nc.contact2)
ncatt_put(nc, varid=0, attname='contact3', attval=nc.contact3)
ncatt_put(nc, varid=0, attname='history', attval=nc.history)

nc_sync(nc)
nc_close(nc)

##******************************************************************************

nc <- nc_open(template.file, write=TRUE)
pr <- ncvar_get(nc, varid='pr')
pr[] <- missing_value
ncvar_put(nc, varid='pr', vals=pr)

tasmax <- ncvar_get(nc, varid='tasmax')
tasmax[] <- missing_value
ncvar_put(nc, varid='tasmax', vals=tasmax)

tasmin <- ncvar_get(nc, varid='tasmin')
tasmin[] <- missing_value
ncvar_put(nc, varid='tasmin', vals=tasmin)

wind <- ncvar_get(nc, varid='wind')
wind[] <- missing_value
ncvar_put(nc, varid='wind', vals=wind)

##******************************************************************************

nc_close(nc)
