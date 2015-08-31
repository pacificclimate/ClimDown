##******************************************************************************
# Bias Corrected Climate Imprint (BCCI) downscaling algorithm
# Alex Cannon (acannon@uvic.ca)
# *All wind related lines of code - commented out by Arelia Werner (wernera@uvic.ca)
##******************************************************************************
# Bilinearly interpolate daily GCM anomalies to fine-scale grid and
# superimpose back onto fine-scale monthly climatologies
##******************************************************************************
##Modified to copy and prepare the new netcdf file

ptm <- proc.time()

library(fields)
library(ncdf4)

code.dir <- '/home/ssobie/stat.downscaling/code/QPQM/BCCI/'

args <- commandArgs(trailingOnly=TRUE)
for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
}

config <- paste(code.dir,'BCCI.set.config',sep='')
print(readLines(config))
source(config)
  
  system(paste('cp ', template.file,'_tasmin_only.nc', ' ', output.dir, 'imp.', output.file,
               output.suffix, '_tasmin.nc', sep=''))
  nc <- nc_open(paste(output.dir, 'imp.', output.file, output.suffix, '_tasmin.nc',
                      sep=''), write=TRUE)
  
  ncatt_put(nc, varid=0, attname='title', attval=output.nc.title)
  ncatt_put(nc, varid=0, attname='institution',
            attval=output.nc.institution)
  ncatt_put(nc, varid=0, attname='source', attval=output.nc.source)
  ncatt_put(nc, varid=0, attname='input_data', attval=output.nc.input_data)
  ncatt_put(nc, varid=0, attname='reference', attval=output.nc.reference)
  ncatt_put(nc, varid=0, attname='project_id', attval=output.nc.project_id)
  ncatt_put(nc, varid=0, attname='experiment_id',
            attval=output.nc.experiment_id)
  ncatt_put(nc, varid=0, attname='version', attval=output.nc.version)
  ncatt_put(nc, varid=0, attname='version_comment',
            attval=output.nc.version_comment)
  ncatt_put(nc, varid=0, attname='contact1', attval=output.nc.contact1)
  ncatt_put(nc, varid=0, attname='contact2', attval=output.nc.contact2)
  ncatt_put(nc, varid=0, attname='contact3', attval=output.nc.contact3)
  ncatt_put(nc, varid=0, attname='history', attval=output.nc.history)
  ncatt_put(nc, varid='time', attname='units', attval=output.nc.time.units)
  ncatt_put(nc, varid='time', attname='calendar', attval=output.nc.time.calendar)
  nc_sync(nc)
  
nc_close(nc)	

print('Elapsed time')
print(proc.time() - ptm)

##******************************************************************************

rm(list=ls())



