################################################################################
# Quantile perturbation quantile mapping for bias correction of ratio and
# interval data. Alex Cannon (acannon@uvic.ca)
################################################################################

library(ncdf4)
library(abind)
library(PCICt)
library(iterators)

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

mk.factor.set <- function(o.c.dates, m.c.dates, m.p.dates,
                          multiyear=FALSE, seasonal=TRUE,
                          n.multiyear=10, expand.multiyear=TRUE) {
    if (seasonal) {
      stop("Seasonal sliding window is not presently implenented")
    }

    if (multiyear) {
        list(
            oc = mk.monthly.factor(o.c.dates),
            mc = mk.monthly.factor(m.c.dates),
            mp = interaction(mk.multiyear.factor(m.p.dates, n.multiyear, expand.multiyear), mk.monthly.factor(m.p.dates), sep='-')
        )
    } else{
        list(
            oc = mk.monthly.factor(o.c.dates),
            mc = mk.monthly.factor(m.c.dates),
            mp = interaction(mk.annual.factor(m.p.dates), mk.monthly.factor(m.p.dates), sep='-')
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
    require(foreach)

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

    return(unsplit(mhat.list, m.p.factor))

    mhat.list <- lapply(years, function(year) {
        cases.p <- mhat.p <- c()
        for(month in months){
            # Current month
            cases.o.c <- which(dates.o.c[,2] %in% month)
            cases.m.c <- which(dates.m.c[,2] %in% month)
            cases.window <- which((dates.m.p[,1] %in% start:end) &
                                  (dates.m.p[,2] %in% month))
            if(seasonal) {
                # Previous month-1
                month.m1 <- ifelse((month-1)==(min(months)-1),
                                   max(months), month-1)
                cases.o.c.m1 <- which(dates.o.c[,2] %in% month.m1)
                cases.m.c.m1 <- which(dates.m.c[,2] %in% month.m1)
                cases.window.m1 <- which((dates.m.p[,1] %in% start:end) &
                                         (dates.m.p[,2] %in% month.m1))
                # Next month+1
                month.p1 <- ifelse((month+1)==(max(months)+1),
                                   min(months), month+1)
                cases.o.c.p1 <- which(dates.o.c[,2] %in% month.p1)
                cases.m.c.p1 <- which(dates.m.c[,2] %in% month.p1)
                cases.window.p1 <- which((dates.m.p[,1] %in% start:end) &
                                         (dates.m.p[,2] %in% month.p1))
                # Bias correction for centre month of seasonal window
                qpqm <- QPQM(o.c=c(o.c[cases.o.c], o.c[cases.o.c.m1],
                             o.c[cases.o.c.p1]), m.c=c(m.c[cases.m.c],
                             m.c[cases.m.c.m1], m.c[cases.m.c.p1]),
                             m.p=c(m.p[cases.window], m.p[cases.window.m1],
                             m.p[cases.window.p1]), ratio=ratio, trace=trace,
                             jitter.factor=jitter.factor, n.tau=n.tau)
                cases.mhat.p <- which(dates.m.p[cases.window,1] %in% year)
                mhat.year <- qpqm$mhat.p[1:length(m.p[cases.window])
                                        ][cases.mhat.p]
            }
            cases.p <- c(cases.p, which((dates.m.p[,1] %in% year) &
                         (dates.m.p[,2] %in% month)))
            mhat.p <- c(mhat.p, mhat.year)
        }
        list(cases.p=cases.p, mhat.p=mhat.p)
    })
}


qpqm.netcdf.wrapper <- function(obs.file, gcm.file, out.file, varname='tasmax') {
    ptm <- proc.time()

    # FIXME: Parametrize all of this
    n.window <- 1
    multiyear <- TRUE
    expand.multiyear <- TRUE
    n.multiyear <- 30
    trace <- 0.005
    jitter.factor <- 0.01
    cstart <- as.POSIXct('1971-01-01', tz='GMT')
    cend <- as.POSIXct('2000-12-31', tz='GMT')
    n.chunks <- 500

    tau <- list(pr=1001, tasmax=101, tasmin=101)

    seasonal <- list(pr=TRUE, tasmax=TRUE, tasmin=TRUE)
    ratio <- list(pr=TRUE, tasmax=FALSE, tasmin=FALSE)

    # Read in the input and output files
    gcm <- nc_open(gcm.file)
    obs <- nc_open(obs.file)

    cat('Creating output file', out.file, '\n')
    dims <- gcm$var[[varname]]$dim
    vars <- ncvar_def(varname, gcm$var[[varname]]$units, dims)
    out <- nc_create(out.file, vars)

    lat <- gcm$dim$lat$vals
    lon <- gcm$dim$lon$vals

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

    gcm.time <- compute.time.stats(gcm, cstart)
    obs.time <- compute.time.stats(obs, cstart, cend)
    # indices for gcm time that are within the observational time range
    gcm.obs.subset.i <- gcm.time$vals > as.PCICt(cstart, attr(gcm.time$vals, 'cal')) & gcm.time$vals < as.PCICt(cend, attr(gcm.time$vals, 'cal'))

    # Calculate the time factors outside of the main spatial loop
    time.factors <- mk.factor.set(obs.time$vals,
                                  gcm.time$vals[gcm.obs.subset.i],
                                  gcm.time$vals,
                                  multiyear=multiyear, seasonal=FALSE,
                                  n.multiyear=n.multiyear, expand.multiyear=expand.multiyear
                                  )

    na.gcm <- rep(NA, gcm.time$n)

    # Calculate out to split up the chunks for reading.
    # We have to read all of time and as many rows as possible
    # The spatial domains are the same
    chunk.size <- optimal.chunk.size((gcm.time$n + obs.time$n) * length(lat))
    chunk.lon.indices <- chunk.indices(length(lon), chunk.size)
    n.chunks <- length(chunk.lon.indices)

    for (chunk in chunk.lon.indices) {
        print(sort( sapply(ls(),function(x){object.size(get(x))})))

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

        m.p.chunk <- mapply(
                       function(i, j) {
                         print(paste(i, ',', j, '/', xn, ',', yn))
                         o.c <- o.c.chunk[i, j, ]
                         m.p <- m.p.chunk[i, j, ]

                         # FIXME: 
                         if(all(is.na(o.c), is.na(m.p))) {
                           na.gcm
                         } else {
                           # consider the modeled values during the observed period separately
                           # m.c <- m.p[gcm.obs.subset.i]

                           tQPQM(o.c=o.c, m.c=m.p[gcm.obs.subset.i], m.p=m.p,
                                 time.factors$oc,
                                 time.factors$mc,
                                 time.factors$mp,
                                 ratio=ratio[[varname]], trace=trace,
                                 jitter.factor=jitter.factor,
                                 n.tau=tau[[varname]])
                         }
                       },
                       ij$i, ij$j)

        dim(m.p.chunk) <- c(gcm.time$n, xn, yn)
        print(dim(m.p.chunk))
        m.p.chunk <- aperm(m.p.chunk, c(2, 3, 1))
        print(dim(m.p.chunk))
        print(paste("Writing chunk", chunk['start'], '-', chunk['stop'], '/', length(lon), 'to file:', out$filename))
        ncvar_put(nc=out, varid=varname, vals=m.p.chunk,
                  start=c(chunk['start'], 1, 1), count=dim(m.p.chunk))
        #print(object.size(x=lapply(ls(), get)), units="Mb")
        rm(o.c.chunk, m.p.chunk)
        gc()
        #print(sort( sapply(ls(),function(x){object.size(get(x))})))
    }


    nc_close(gcm)
    nc_close(obs)
    nc_close(out)

    print('Elapsed Time')
    print(proc.time() - ptm)
}
################################################################################
