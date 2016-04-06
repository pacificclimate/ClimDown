################################################################################
# Quantile perturbation quantile mapping for bias correction of ratio and
# interval data. Alex Cannon (acannon@uvic.ca)
################################################################################

QPQM <- function(o.c, m.c, m.p, ratio=TRUE, trace=0.05, jitter.factor=0.01,
                 n.tau=NULL)
{
    # Quantile perturbation quantile mapping bias correction:
    # o = vector of observed values; m = vector of modelled values
    # c = current period;  p = projected period
    # ratio = TRUE --> preserve relative trends in a ratio variable
    # trace = 0.05 --> treat values below trace as left censored
    # jitter.factor = 0.01 --> jittering to accomodate ties
    # n.tau = NULL --> number of empirical quantiles (NULL=sample length)
    #
    # tau.m-p = F.m-p(x.m-p)
    # delta.m = x.m-p {/,-} F.m-c^-1(tau.m-p)
    # xhat.m-p = F.o-c^-1(tau.m-p) {*,+} delta.m
    #
    # Apply a small amount of jitter to accomodate ties due to limited
    # measurement precision
    if (all(is.na(m.p))) {
      return(list(mhat.c=NULL, mhat.p=rep(NA, length(m.p))))
    }
    o.c <- jitter(o.c, jitter.factor)
    m.c <- jitter(m.c, jitter.factor)
    m.p <- jitter(m.p, jitter.factor)
    # For ratio data, treat exact zeros as left censored values less than trace
    if(ratio){
        epsilon <- .Machine$double.eps
        o.c[o.c < trace & !is.na(o.c)] <- runif(sum(o.c < trace, na.rm=TRUE), min=epsilon, max=trace)
        m.c[m.c < trace & !is.na(m.c)] <- runif(sum(m.c < trace, na.rm=TRUE), min=epsilon, max=trace)
        m.p[m.p < trace & !is.na(m.p)] <- runif(sum(m.p < trace, na.rm=TRUE), min=epsilon, max=trace)
      }
    # Calculate empirical quantiles using Weibull plotting position
    n <- max(length(o.c), length(m.c), length(m.p))
    if(is.null(n.tau)) n.tau <- n
    tau <- seq(1/(n+1), n/(n+1), length=n.tau)
    quant.o.c <- quantile(o.c, tau, type=6, na.rm=TRUE)
    quant.m.c <- quantile(m.c, tau, type=6, na.rm=TRUE)
    quant.m.p <- quantile(m.p, tau, type=6, na.rm=TRUE)
    # Apply QPQM bias correction
    tau.m.p <- approx(quant.m.p, tau, m.p, rule=2)$y    
    if(ratio){
        delta.m <- m.p/approx(tau, quant.m.c, tau.m.p, rule=2)$y
        mhat.p <- approx(tau, quant.o.c, tau.m.p, rule=2)$y*delta.m
    } else{
        delta.m <- m.p - approx(tau, quant.m.c, tau.m.p, rule=2)$y
        mhat.p <- approx(tau, quant.o.c, tau.m.p, rule=2)$y + delta.m
    }
    mhat.c <- approx(quant.m.c, quant.o.c, m.c, rule=2)$y
    # For ratio data, set values less than trace to zero
    if(ratio){
        mhat.c[mhat.c < trace] <- 0
        mhat.p[mhat.p < trace] <- 0
    }
    list(mhat.c=mhat.c, mhat.p=mhat.p)
}

mk.multiyear.factor <- function(dates, block.size, expand.multiyear=TRUE) {
    years <- as.numeric(format(dates, '%Y'))

    block.factor <- factor(years - years %% block.size)

    if(expand.multiyear) {
        # Fold incomplete data from the final block into the previous
        # complete multi-year block
        multiyear.lengths <- tapply(block.factor, block.factor, length)

        n.thresh <- max(multiyear.lengths[-length(multiyear.lengths)])/2
        # Is the final block less than half the length of other blocks?
        if (multiyear.lengths[nlevels(block.factor)] < n.thresh) {
            # Assign the final level to the value of the penultimate level
            levels(block.factor)[nlevels(block.factor)] <- levels(block.factor)[nlevels(block.factor)-1]
        }
    }

    block.factor
}

mk.annual.factor <- function(dates) {
    factor(format(dates, '%Y'))
}

mk.monthly.factor <- function(dates) {
  factor(format(dates, '%m'))
}

mk.seasonal.factor <- function(dates) {
    cal <- attr(dates, 'cal')
    # PCICt can possibly omit the "_day" on the cal attribute. Add it back.
    cal <- sub('^(360|365)$', '\\1_day', cal)
    if (cal == "proleptic_gregorian") { cal <- NULL}
    jday <- as.integer(format(dates, '%j'))
    years <- as.integer(format(dates, '%Y'))
    mkseas(jday, 'DJF', calendar=cal, year=years)
}

mk.factor.set <- function(o.c.dates, m.c.dates, m.p.dates,
                          multiyear=FALSE, seasonal=TRUE,
                          n.multiyear=10, expand.multiyear=TRUE) {
    if (seasonal) {
        mk.factor <- mk.seasonal.factor
    } else {
        mk.factor <- mk.monthly.factor
    }

    if (multiyear) {
        list(
            oc = mk.factor(o.c.dates),
            mc = mk.factor(m.c.dates),
            mp = interaction(mk.multiyear.factor(m.p.dates, n.multiyear, expand.multiyear), mk.factor(m.p.dates), sep='-')
        )
    } else{
        list(
            oc = mk.factor(o.c.dates),
            mc = mk.factor(m.c.dates),
            mp = interaction(mk.annual.factor(m.p.dates), mk.factor(m.p.dates), sep='-')
        )
    }
}

tQPQM <- function(o.c, m.c, m.p,
                  o.c.factor, m.c.factor, m.p.factor,
                  n.window=30, ratio=TRUE, trace=0.05,
                  jitter.factor=0.01, n.tau=NULL)
{
    # Apply QPQM bias correction over 3-month moving windows (seasonal=TRUE)
    # or single months, both over (sliding) blocks of n.window years or
    # multi-year chunks of years
    # o = vector of observed values; m = vector of modelled values
    # c = current period;  p = projected period
    # *.factor = date factors for grouping values together, levels should
    # be the same in each of the three factors
    # ratio = TRUE --> preserve relative trends in a ratio variable
    # trace = 0.05 --> treat values below trace as left censored
    # jitter.factor = 0.01 --> jittering to accomodate ties
    # seasonal = TRUE --> apply over sliding 3-month windows
    # multiyear = FALSE --> apply over multi-year rather than 1-yr chunks
    # n.multiyear = 10 --> if multiyear==TRUE, apply to n.multiyear chunks
    # expand.multiyear --> fold incomplete multi-year block into previous
    # n.tau = NULL --> number of empirical quantiles (NULL=sample length)

    # We'll run the computation on blocks of data that consist of the same month/season
    # and maybe across multiple (e.g. 30) years (if multiyear is TRUE)
    sections <- split(m.p, m.p.factor)
    months <- sub('^[0-9]+-', '', names(sections))
    oc.by.month <- split(o.c, o.c.factor)
    mc.by.month <- split(m.c, m.c.factor)
    mhat.list <- mapply(
                   function(mp.subset, month) {
                     qpqm <- QPQM(o.c=oc.by.month[[month]], m.c=mc.by.month[[month]],
                                  m.p=mp.subset, ratio=ratio,
                                  trace=trace, jitter.factor=jitter.factor,
                                  n.tau=n.tau)
                     qpqm$mhat.p
                   },
                   mp.subset=sections,
                   months
                   )

    unsplit(mhat.list, m.p.factor)
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

#' @title High-level wrapper for Quantile perturbation quantile mapping (QPQM)
#'
#' @description This function performs the QPQM algorithm on a
#' cell-by-cell basis for each cell in the spatial domain of the
#' inputted high-res gridded observations. It uses the gridded
#' observations plus the GCM-based output of BCCI as input to the
#' algorithm and then performs a quantile perturbation/quantile
#' mapping bias correction. The output is written out to out.file.
#' 
#' @param obs.file Filename of high-res gridded historical observations
#' @param gcm.file Filename of GCM simulations interpolated to the obs.file grid
#' @param out.file Filename to create (or overwrite) with the bias corrected outputs
#' @param varname Name of the NetCDF variable to downscale (e.g. 'tasmax')
#' @return NULL
#'
#' @export
qpqm.netcdf.wrapper <- function(obs.file, gcm.file, out.file, varname='tasmax') {
    ptm <- proc.time()

    cstart <- getOption('cstart')
    cend <- getOption('cend')

    print("Opening the input files and reading metadata")
    gcm <- nc_open(gcm.file)
    obs <- nc_open(obs.file)

    lat <- gcm$dim$lat$vals
    lon <- gcm$dim$lon$vals

    gcm.time <- compute.time.stats(gcm, cstart)
    obs.time <- compute.time.stats(obs, cstart, cend)
    # indices for gcm time that are within the observational time range
    gcm.obs.subset.i <- gcm.time$vals > as.PCICt(cstart, attr(gcm.time$vals, 'cal')) & gcm.time$vals < as.PCICt(cend, attr(gcm.time$vals, 'cal'))

    print("Calculating the time factors outside in order to subdivide the problem space")
    time.factors <- mk.factor.set(obs.time$vals,
                                  gcm.time$vals[gcm.obs.subset.i],
                                  gcm.time$vals,
                                  multiyear=getOption('multiyear'),
                                  seasonal=getOption('seasonal')[[varname]],
                                  n.multiyear=getOption('multiyear.window.length'),
                                  expand.multiyear=getOption('expand.multiyear')
                                  )

    cat('Creating output file', out.file, '\n')
    # FIXME: The GCM time needs to be clipped to cstart
    dims <- gcm$var[[varname]]$dim
    vars <- ncvar_def(varname, gcm$var[[varname]]$units, dims)
    out <- nc_create(out.file, vars)

    na.gcm <- rep(NA, gcm.time$n)

    # Calculate out to split up the chunks for reading.
    # We have to read all of time and as many rows as possible
    # The spatial domains are the same
    chunk.size <- optimal.chunk.size((gcm.time$n + obs.time$n) * length(lat))
    chunk.lon.indices <- chunk.indices(length(lon), chunk.size)
    n.chunks <- length(chunk.lon.indices)

    be <- start.par.backend()
    # I/O loop (read, compute, write)
    for (chunk in chunk.lon.indices) {

        print(paste('Bias correcting', varname, 'longitudes', chunk['start'], '-', chunk['stop'], '/', length(lon)))
        print(paste("Reading longitudes", chunk['start'], '-', chunk['stop'], '/', length(lon), 'from file:', obs$filename))
        print(paste("... and reading latitudes 1 -", length(lat), '/', length(lat)))

        o.c.chunk <- ncvar_get(obs, start=c(chunk['start'], 1, obs.time$t0),
                               count=c(chunk['length'], -1, obs.time$n),
                               varid=varname, collapse_degen=FALSE)

        print(paste("Reading longitudes", chunk['start'], '-', chunk['stop'], '/', length(lon), 'from file:', gcm$filename))
        print(paste("... and reading latitudes 1 -", length(lat), '/', length(lat)))

        m.p.chunk <- ncvar_get(gcm, start=c(chunk['start'], 1, gcm.time$t0),
                               count=c(chunk['length'], -1, gcm.time$n),
                               varid=varname, collapse_degen=FALSE)

        xn <- dim(o.c.chunk)[1]
        yn <- dim(o.c.chunk)[2]
        ij <- expand.grid(i=seq(xn), j=seq(yn))

        ncells <- xn*yn

        # Compute loop. Cell major. Do something to a timeseries for each cell.
        print(paste("Computing QPQM on", ncells, "cells"))
        m.p.chunk <- foreach(o.c=split(o.c.chunk, 1:ncells),
                             m.p=split(m.p.chunk, 1:ncells),
                             .combine=c, .multicombine=TRUE, .inorder=TRUE,
                             .export=c('na.gcm', 'tQPQM', 'varname', 'gcm.obs.subset.i', 'time.factors', 'QPQM')) %loop% {

          if(all(is.na(o.c), is.na(m.p))) {
            na.gcm
          } else {
            # consider the modeled values during the observed period separately
            # m.c <- m.p[gcm.obs.subset.i]
            tQPQM(o.c=o.c, m.c=m.p[gcm.obs.subset.i], m.p=m.p,
                  time.factors$oc,
                  time.factors$mc,
                  time.factors$mp,
                  ratio=getOption('ratio')[[varname]],
                  trace=getOption('trace'),
                  jitter.factor=getOption('jitter.factor'),
                  n.tau=getOption('tau')[[varname]])
          }
        }

        dim(m.p.chunk) <- c(gcm.time$n, xn, yn)
        m.p.chunk <- aperm(m.p.chunk, c(2, 3, 1))

        print(paste("Writing chunk", chunk['start'], '-', chunk['stop'], '/', length(lon), 'to file:', out$filename))
        ncvar_put(nc=out, varid=varname, vals=m.p.chunk,
                  start=c(chunk['start'], 1, 1), count=dim(m.p.chunk))
        rm(o.c.chunk, m.p.chunk)
        gc()
    }
    stop.par.backend(be)

    print("Closing up the input and output files")
    nc_close(gcm)
    nc_close(obs)
    nc_close(out)

    print('Elapsed Time')
    print(proc.time() - ptm)
}
################################################################################
