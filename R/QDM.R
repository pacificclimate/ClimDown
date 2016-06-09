reorder <- function(x,ix) {
  rx <- x[ix]
  return(rx)
}

#' @title High-level NetCDF wrapper for Quantile Delta Mapping (QDM)
#'
#' @description All files (save for the analogues_file) should have the same spatial domain.
#'
#' @param qpqm.file The output file from the QPQM script
#' @param obs.file Filename of high-res gridded historical observations
#' @param analogues Temporal analogues... describe this more
#' @param out.file The file to create (or overwrite) with the final NetCDF output
#' @param varname Name of the NetCDF variable to downscale (e.g. 'tasmax')
#'
#' @export
qdm.netcdf.wrapper <- function(qpqm.file, obs.file, analogues, out.file, varname='tasmax') {
    ptm <- proc.time()

    qpqm.nc <- nc_open(qpqm.file)
    obs.nc <- nc_open(obs.file)


    cat('Creating output file', out.file, '\n')
    dims <- qpqm.nc$var[[varname]]$dim
    vars <- ncvar_def(varname, getOption('target.units')[varname], dims)
    out.nc <- nc_create(out.file, vars)

    nlat <- qpqm.nc$dim$lat$len
    nlon <- qpqm.nc$dim$lon$len

    qpqm.time <- compute.time.stats(qpqm.nc)

    n.qpqm <- qpqm.nc$dim$time$len
    nt <- qpqm.time$n

    ncells <- nlat * nlon
    chunk.size <- optimal.chunk.size(ncells)
    fact <- chunk.month.factor(qpqm.time$vals, chunk.size)

    # Do I/O in large chunks
    indices <- lapply(
        split(1:qpqm.time$n, fact),
        function(x) {c('start'=min(x), 'stop'=max(x), 'length'=length(x))}
    )

    for (index in indices) {
        i_0 <- index['start']
        i_n <- index['stop']
        ni <- index['length']

        print(paste("Processing steps", i_0, "-", i_n, "/", nt))

        # Read all space for the time chunk
        var.qpqm <- ncvar_get(qpqm.nc, varid=varname, start=c(1,1,i_0), count=c(nlon,nlat,ni))

        date.sub <- qpqm.time$vals[i_0:i_n]
        month.factor <- as.factor(format(date.sub, '%Y-%m'))

        print(paste("Applying analogues to timesteps", i_0, "-", i_n, "/", nt))
        var.bcca <- foreach(
            ti=analogues$indices[i_0:i_n],
            wi=analogues$weights[i_0:i_n],
            .export=c('obs.nc', 'varname'),
            .final=function(x) {
                array(unlist(x), dim=c(ncells, nt))
            }
            ) %dopar% {
                apply.analogues.netcdf(ti, wi, obs.nc, varname)
            }

        dqm <- foreach(
            bcca=split(var.bcca, rep(month.factor, each=ncells)),
            qpqm=split(var.qpqm, rep(month.factor, each=ncells)),
            .multicombine=TRUE,
            .inorder=TRUE,
            .final=function(x) {
                array(
                    unsplit(x, rep(month.factor, each=ncells)),
                    dim=c(nlon, nlat, ni)
                )
            }) %dopar% {
                bcca <- jitter(bcca, 0.01)
                ndays <- length(bcca) / ncells
                dim(bcca) <- c(nlon, nlat, ndays)
                ranks <- apply(bcca, 1:2, rank, ties.method='average')
                ## Reorder the days of qpqm
                array(
                    mapply(
                        reorder,
                        split(qpqm, rep(1:ndays, each=ncells)),
                        split(ranks, 1:ndays)
                    ),
                    dim=c(nlon, nlat, ndays)
                )
            }
        print(paste("Writing steps", i_0, "-", i_n, "/", nt, "to file", out.file))
        ncvar_put(nc=out.nc, varid=varname, vals=dqm,
                  start=c(1, 1, i_0), count=c(-1, -1, ni))
    }

    nc_close(qpqm.nc)
    nc_close(obs.nc)
    nc_close(out.nc)

    print('Elapsed time')
    print(proc.time() - ptm)
}
