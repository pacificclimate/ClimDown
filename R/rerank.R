reorder <- function(x,ix) {
  rx <- sort(x)[ix]
  return(rx)
}

#' @title High-level NetCDF wrapper for Quantile Reranking
#'
#' @description All files (save for the analogues_file) should have the same spatial domain.
#'
#' @param qdm.file The output file from the QDM script
#' @param obs.file Filename of high-res gridded historical observations
#' @param analogues Temporal analogues... describe this more
#' @param out.file The file to create (or overwrite) with the final NetCDF output
#' @param varname Name of the NetCDF variable to downscale (e.g. 'tasmax')
#'
#' @references Schefzik, R., Thorarinsdottir, T. L., & Gneiting, T. (2013). Uncertainty quantification in complex simulation models using ensemble copula coupling. Statistical Science, 28(4), 616-640.
#'
#' Wilks, D. S. (2015). Multivariate ensemble Model Output Statistics using empirical copulas. Quarterly Journal of the Royal Meteorological Society, 141(688), 945-952.
#' @export
rerank.netcdf.wrapper <- function(qdm.file, obs.file, analogues, out.file, varname='tasmax') {
    ptm <- proc.time()

    qdm.nc <- nc_open(qdm.file)
    obs.nc <- nc_open(obs.file)


    cat('Creating output file', out.file, '\n')
    dims <- qdm.nc$var[[varname]]$dim
    vars <- ncvar_def(varname, getOption('target.units')[varname], dims)
    out.nc <- nc_create(out.file, vars)

    nlat <- qdm.nc$dim$lat$len
    nlon <- qdm.nc$dim$lon$len

    qdm.time <- compute.time.stats(qdm.nc)

    n.qdm <- qdm.nc$dim$time$len
    nt <- qdm.time$n

    ncells <- nlat * nlon
    chunk.size <- optimal.chunk.size(ncells)
    fact <- chunk.month.factor(qdm.time$vals, chunk.size)

    # Do I/O in large chunks
    indices <- lapply(
        split(1:qdm.time$n, fact),
        function(x) {c('start'=min(x), 'stop'=max(x), 'length'=length(x))}
    )

    for (index in indices) {
        i_0 <- index['start']
        i_n <- index['stop']
        ni <- index['length']

        print(paste("Processing steps", i_0, "-", i_n, "/", nt))

        # Read all space for the time chunk
        var.qdm <- ncvar_get(qdm.nc, varid=varname, start=c(1,1,i_0), count=c(nlon,nlat,ni))

        date.sub <- qdm.time$vals[i_0:i_n]
        month.factor <- as.factor(format(date.sub, '%Y-%m'))

        print(paste("Applying analogues to timesteps", i_0, "-", i_n, "/", nt))
        var.ca <- mapply(
            function(ti, wi) {
                apply.analogues.netcdf(ti, wi, obs.nc, varname)
            },
            analogues$indices[i_0:i_n],
            analogues$weights[i_0:i_n]
        )
        var.ca <- positive_pr(var.ca, varname)
        by.month <- rep(month.factor, each=ncells)

        dqm <- foreach(
            ca=split(var.ca, by.month),
            qdm=split(var.qdm, by.month),
            .multicombine=TRUE,
            .inorder=TRUE,
            .final=function(x) {
                array(
                    unsplit(x, by.month),
                    dim=c(nlon, nlat, ni)
                )
            }) %dopar% {
                ca <- jitter(ca, 0.01)
                ndays <- length(ca) / ncells
                by.cell <- rep(1:ncells, times=ndays)
                ca <- split(ca, by.cell)
                ranks <- lapply(ca, rank, ties.method='average')
                ## Reorder the days of qdm on cell at a time
                array(
                    unsplit(
                        mapply(
                            reorder,
                            split(qdm, by.cell),
                            ranks,
                            SIMPLIFY=F
                        ),
                        by.cell
                    ),
                    dim=c(nlon, nlat, ndays)
                )
            }
        print(paste("Writing steps", i_0, "-", i_n, "/", nt, "to file", out.file))
        ncvar_put(nc=out.nc, varid=varname, vals=dqm,
                  start=c(1, 1, i_0), count=c(-1, -1, ni))
        rm(var.qdm, var.ca, dqm)
        gc()
    }

    nc_close(qdm.nc)
    nc_close(obs.nc)
    nc_close(out.nc)

    print('Elapsed time')
    print(proc.time() - ptm)
}
