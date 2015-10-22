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
    o.c <- jitter(o.c, jitter.factor)
    m.c <- jitter(m.c, jitter.factor)
    m.p <- jitter(m.p, jitter.factor)
    # For ratio data, treat exact zeros as left censored values less than trace
    if(ratio){
        epsilon <- .Machine$double.eps
        o.c[o.c < trace] <- runif(sum(o.c < trace), min=epsilon, max=trace)
        m.c[m.c < trace] <- runif(sum(m.c < trace), min=epsilon, max=trace)
        m.p[m.p < trace] <- runif(sum(m.p < trace), min=epsilon, max=trace)
    }
    # Calculate empirical quantiles using Weibull plotting position
    n <- max(length(o.c), length(m.c), length(m.p))
    if(is.null(n.tau)) n.tau <- n
    tau <- seq(1/(n+1), n/(n+1), length=n.tau)
    quant.o.c <- quantile(o.c, tau, type=6)
    quant.m.c <- quantile(m.c, tau, type=6)
    quant.m.p <- quantile(m.p, tau, type=6)
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

tQPQM <- function(o.c, m.c, m.p, dates.o.c, dates.m.c, dates.m.p,
                  n.window=30, ratio=TRUE, trace=0.05,
                  jitter.factor=0.01, seasonal=TRUE, multiyear=FALSE,
                  n.multiyear=10, expand.multiyear=TRUE, n.tau=NULL)
{
    # Apply QPQM bias correction over 3-month moving windows (seasonal=TRUE)
    # or single months, both over (sliding) blocks of n.window years or
    # multi-year chunks of years
    # o = vector of observed values; m = vector of modelled values
    # dates = matrix of dates with integer year, month, day columns
    # c = current period;  p = projected period
    # n.window = window length for blocks
    # ratio = TRUE --> preserve relative trends in a ratio variable
    # trace = 0.05 --> treat values below trace as left censored
    # jitter.factor = 0.01 --> jittering to accomodate ties
    # seasonal = TRUE --> apply over sliding 3-month windows
    # multiyear = FALSE --> apply over multi-year rather than 1-yr chunks
    # n.multiyear = 10 --> if multiyear==TRUE, apply to n.multiyear chunks
    # expand.multiyear --> fold incomplete multi-year block into previous
    # n.tau = NULL --> number of empirical quantiles (NULL=sample length)
    require(foreach)
    months <- unique(dates.m.p[,2])
    print(paste('i is: ',i,sep=''))
    print(paste('j is: ',j,sep=''))
    if(multiyear){
        # Apply QPQM to multi-year blocks of length n.multiyear
        dates.m.p[,1] <- dates.m.p[,1]-(dates.m.p[,1]%%n.multiyear)
        if(expand.multiyear){
            # Fold incomplete data from the final block into the previous
            # complete multi-year block
            multiyear.lengths <- sapply(split(dates.m.p[,1],
                                        dates.m.p[,1] -
                                        (dates.m.p[,1]%%n.multiyear)),
                                        length)
            n.thresh <- max(multiyear.lengths[-length(multiyear.lengths)])/2
            print(n.thresh)
            incomplete <- names(which(multiyear.lengths < n.thresh))
            if(length(incomplete) > 0){
                dates.m.p[dates.m.p[,1]==incomplete,1] <-
                    max(dates.m.p[dates.m.p[,1]!=incomplete,1])
            }
        }
        chunks <- unique(dates.m.p[,1])
        years <- seq_along(chunks)
        for(i in seq_along(years))
            dates.m.p[dates.m.p[,1]==chunks[i],1] <- i
    } else{
        years <- unique(dates.m.p[,1])
    }
    dplus <- dminus <- floor(n.window/2)
    dplus <- ifelse((dplus+dminus)==n.window, dplus-1, dplus)

    mhat.list <- lapply(years, function(year) {
        start <- max(years[1], year-dminus)
        end <- min(years[length(years)], year+dplus)
        if((end-start+1) < n.window){
            dyears <- (n.window-(end-start+1))
            if(start==years[1]){
                end <- end + dyears
            } else{
                start <- start - dyears
            }
        }
        cases.p <- mhat.p <- c()
        for(month in months){
            # Current month
            cases.o.c <- which(dates.o.c[,2] %in% month)
            cases.m.c <- which(dates.m.c[,2] %in% month)
            cases.window <- which((dates.m.p[,1] %in% start:end) &
                                  (dates.m.p[,2] %in% month))
            if(seasonal){
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
            } else{
                # Bias correction by individual month
                qpqm <- QPQM(o.c=o.c[cases.o.c], m.c=m.c[cases.m.c],
                             m.p=m.p[cases.window], ratio=ratio,
                             trace=trace, jitter.factor=jitter.factor,
                             n.tau=n.tau)
                cases.mhat.p <- which(dates.m.p[cases.window,1] %in% year)
                mhat.year <- qpqm$mhat.p[cases.mhat.p]
            }
            cases.p <- c(cases.p, which((dates.m.p[,1] %in% year) &
                         (dates.m.p[,2] %in% month)))
            mhat.p <- c(mhat.p, mhat.year)
        }
        list(cases.p=cases.p, mhat.p=mhat.p)
    })
    mhat.p <- m.p*NA
    for(i in seq_along(mhat.list))
        mhat.p[mhat.list[[i]]$cases.p] <- mhat.list[[i]]$mhat.p
    mhat.p
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

    pr.n.tau <- 1001
    tasmax.n.tau <- 101
    tasmin.n.tau <- 101

    seasonal <- list(pr=TRUE, tasmax=TRUE, tasmin=TRUE)
    ratio <- seasonal

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

        o.c.chunk <- ncvar_get(obs, start=c(chunk['start'], 1, obs.time$t0),
                               count=c(-1, chunk['length'], obs.time$n),
                               varid=varname, collapse_degen=FALSE)

        print(paste("Reading longitudes", chunk['start'], '-', chunk['stop'], '/', length(lon), 'from file:', gcm$filename))

        m.p.chunk <- ncvar_get(gcm, start=c(chunk['start'], 1, gcm.time$t0),
                               count=c(-1, chunk['length'], gcm.time$n),
                               varid=varname, collapse_degen=FALSE)

        xn <- dim(o.c.chunk)[1]
        yn <- dim(o.c.chunk)[2]
        ij <- t(expand.grid(i=seq(xn), j=seq(yn)))

        m.p.chunk <- foreach(ij=ij,
                             o.c=lapply(ij, function(ij) {o.c.chunk[ij['i'], ij['j'],] } ),
                             m.p=lapply(ij, function(ij) {m.p.chunk[ij['i'], ij['j'],] } ),
                             .inorder=TRUE,
                             .combine=cbind,
                             .multicombine=TRUE
                             #.options.mpi=mpi.options) %do% {
                             ) %do% {

                               # FIXME: 
                               if(all(is.na(o.c), is.na(m.p))) {
                                 na.gcm
                               } else {
                                 # consider the modeled values during the observed period separately
                                 dates.m.c <- gcm.time$vals[gcm.time$vals < max(obs.time$vals)]
                                 m.c <- m.p[gcm.time$vals < max(obs.time$vals)]
                                 tQPQM(o.c=o.c, m.c=m.c, m.p=m.p, dates.o.c=obs.time$vals,
                                       dates.m.c=dates.m.c, dates.m.p=gcm.time$vals,
                                       n.window=n.window, ratio=ratio[[varname]], trace=trace,
                                       jitter.factor=jitter.factor, seasonal=seasonal[[varname]],
                                       multiyear=multiyear, n.multiyear=n.multiyear,
                                       expand.multiyear=expand.multiyear, n.tau=pr.n.tau)
                               }
                             }

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
