##******************************************************************************
# Bias Corrected Constructed Analogue (BCCA) downscaling algorithm
# Based loosely off of code by Alex Cannon <acannon@uvic.ca>
# Rewritten by James Hiebert <hiebert@uvic.ca>

library(ncdf4)
code.dir <- Sys.getenv('CODE_DIR')
source(paste(code.dir, 'bisect.R', sep='/'))
source(paste(code.dir, 'netcdf.calendar.R', sep='/'))

# "--args gcm.file='${gcm.file}' obs.file='${obs.file}' varid='${varid}'"
args <- commandArgs(trailingOnly=TRUE)
for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
}

# Input is a factor of which maps fine-scale obs cells to large-scale gcm cells
# The factor should be of length x * y
# and a 3d array of obs (x, y, time)
aggregate.obs <- function(cell.factor, obs) {
  apply(obs, 3, function(x) tapply(x, cell.factor, mean, trim=0.1, na.rm=TRUE))
}

# Input cell indicies mapping obs grid to GCM grid
# and a 3d array of obs (x, y, time)
# Output is 3d array, gcmx x gcm y x time
aggregate.obs.to.gcm.grid <- function(xi, yi, xn, yn, obs) {
  cell.number <- xi * max(yi) + yi
  cell.factor <- factor(cell.number, unique(as.vector(cell.number)))
  # apply preserving time (dim 3)
  # so for each time step, aggregate (take the mean) according to the cell map
  rv <- aggregate.obs(cell.factor, obs)
  ti <- dim(obs)[3]
  dim(rv) <- c(xn, yn, ti)
  return(rv)
}

# Takes a vector length and chunk size
# returns a list of (start, stop, length)
chunk.indices <- function(total.size, chunk.size) {
  lapply(
    split(1:total.size, ceiling(1:total.size / chunk.size)),
    function(x) {c('start'=min(x), 'stop'=max(x), 'length'=length(x))}
    )
}

optimal.chuck.size <- function(n.elements, max.GB=10) {
  # 8 byte numerics
  floor(max.GB * 2 ** 30 / 8 / n.elements)
}

##******************************************************************************
# Read fine-scale grid and spatially aggregate to GCM grid
##******************************************************************************
create.aggregates <- function(obs.file, gcm.file, varid) {

  # Read fine-scale and GCM grid dimensions
  nc.obs <- nc_open(obs.file)
  nc.gcm <- nc_open(gcm.file)
  obs.lons <- ncvar_get(nc.obs, 'lon')
  obs.lats <- ncvar_get(nc.obs, 'lat')
  gcm.lons <- ncvar_get(nc.gcm, 'lon')-360
  gcm.lats <- ncvar_get(nc.gcm, 'lat')
  obs.time <- netcdf.calendar(nc.obs)

  # Figure out which GCM grid boxes are associated with each fine-scale grid point
  grid.mapping <- regrid.coarse.to.fine(gcm.lats, gcm.lons, obs.lats, obs.lons)
  xi <- grid.mapping$xi
  yi <- grid.mapping$yi

  xn <- length(unique(as.vector(xi)))
  yn <- length(unique(as.vector(yi)))
  aggregates <- array(dim=c(length(gcm.lons), length(gcm.lats), length(obs.time)))

  chunk.size <- optimal.chuck.size(length(obs.lons) * length(obs.lats))

  chunks <- chunk.indices(nrow(obs.time), chunk.size)
  # Loop over chunks fo time
  for (i in chunks) {
    cat(i['start'], i['stop'], '\n')
    obs <- ncvar_get(nc.obs, varid=varid, start=c(1, 1, i['start']), # get obs for one chunk
                     count=c(-1, -1, i['length']))
    agg <- aggregate.obs.to.gcm.grid(xi, yi, xn, yn, obs)
    aggregates[min(xi):max(xi), min(yi):max(yi), i['start']:i['stop']] <- agg
    rm(obs)
    gc()
  }
  aggregates
}

print('Starting spatial aggregation')
ptm <- proc.time()
aggregates <- create.aggregates(obs.file, gcm.file, varid)
save(aggregates, file=paste(paste(varid, 'aggregate', 'RData', sep='.')))

print('Elapsed time')
print(proc.time() - ptm)
