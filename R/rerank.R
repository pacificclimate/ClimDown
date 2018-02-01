reorder <- function(x,ix) {
  rx <- sort(x)[ix]
  return(rx)
}

## Package check tools can't detect foreach's use of non-standard evaluation
## Ensure that they skip these variable names
utils::globalVariables(c('ca', 'qdm'))

#' @title High-level NetCDF wrapper for Quantile Reranking
#'
#' @description Quantile Reranking is the final, critical step in the
#'     BCCAQ pipeline. Its purpose is this: since Climate Analogues
#'     (CA) gets its high resolution information by using a linear
#'     combination of historical daily time series for the domain as a
#'     whole, it ends up reintroducing some bias. This is because the
#'     quantile mapping bias correction step was performed only at
#'     course resolution (of the GCM). Quantile Reranking fixes this
#'     by re-applying a simple quantile mapping bias correction at
#'     each grid box. The advantage of doing this as a final step is
#'     that the downscaling method retains the primary advantage of
#'     BCCA: high spatial consistency (e.g. when a storm or a heat
#'     wave hits a specific area, it probably also hits neighboring
#'     areas, etc.).
#'
#' @details All input files (save for the analogues_file) should have
#'     the same spatial domain.
#'
#' @param qdm.file The output file from the QDM script
#' @param obs.file Filename of high-res gridded historical observations
#' @param analogues Temporal analogues. This is a list of two arrays:
#'     the index values and the fractional weights. Each array is the
#'     length of the number of timesteps by k (typically 30).
#' @param out.file The file to create (or overwrite) with the final NetCDF output
#' @param varname Name of the NetCDF variable to downscale (e.g. 'tasmax')
#'
#' @examples
#' \dontrun{
#' options(
#'     calibration.end=as.POSIXct('1972-12-31', tz='GMT')
#' )
#' ci.file <- tempfile(fileext='.nc')
#' ClimDown::ci.netcdf.wrapper('./tiny_gcm.nc', './tiny_obs.nc', ci.file, 'tasmax')
#' qdm.file <<- tempfile(fileext='.nc')
#' ClimDown::qdm.netcdf.wrapper('./tiny_obs.nc', ci.file, qdm.file, 'tasmax')
#' unlink(ci.file)
#' analogues <<- ClimDown::ca.netcdf.wrapper('./tiny_gcm.nc', './tiny_obs.nc')
#' out.file <- tempfile(fileext='.nc')
#' ClimDown::rerank.netcdf.wrapper(qdm.file, './tiny_obs.nc', analogues, out.file, varname='tasmax')
#' unlink(qdm.file)
#' unlink(out.file)
#' }
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
    vars <- ncvar_def(varname, getOption('target.units')[varname], dims, NA)
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
