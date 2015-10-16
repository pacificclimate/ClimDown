################################################################################
# Quantile perturbation quantile mapping for bias correction of ratio and
# interval data. Alex Cannon (acannon@uvic.ca)
################################################################################

library(ncdf4)
library(abind)
library(PCICt)

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

    n.window <- 1
    multiyear <- TRUE
    expand.multiyear <- TRUE
    n.multiyear <- 30
    trace <- 0.005
    jitter.factor <- 0.01
    cstart <- 1971 
    cend <- 2000
    n.chunks <- 500

    pr.n.tau <- 1001
    tasmax.n.tau <- 101
    tasmin.n.tau <- 101

    seasonal <- list(pr=TRUE, tasmax=TRUE, tasmin=TRUE)
    ratio <- seasonal

    # Read in the input and output files
    gcm <- nc_open(gcm.file)
    obs <- nc_open(obs.file)

    cat('Creating input file', out.file, '\n')
    out <- nc_create(out.file, obs$var[[varname]])

    lat <- gcm$dim$lat$vals
    lon <- gcm$dim$lon$vals

    dates.gcm <- as.PCICt(strsplit(gcm$dim$time$units, ' ')[[1]][3],
                          gcm$dim$time$calendar) + gcm$dim$time$vals*86400
    dates.gcm <-  apply(do.call(rbind, strsplit(format(dates.gcm, '%Y %m %d'),
                                                ' ')), 2, as.integer)

    dates.obs <- as.PCICt(strsplit(obs$dim$time$units, ' ')[[1]][3],
                          obs$dim$time$calendar) + obs$dim$time$vals*86400
    dates.obs <-  apply(do.call(rbind, strsplit(format(dates.obs, '%Y %m %d'),
                                                ' ')), 2, as.integer)

    cal.gcm <- c(min(which(dates.gcm[,1]==cstart)), max(which(dates.gcm[,1]==cend)))
    cal.obs <- c(min(which(dates.obs[,1]==cstart)), max(which(dates.obs[,1]==cend)))

    n.gcm <- nrow(dates.gcm)
    n.obs <- diff(cal.obs)+1

    dates.o.c <- dates.obs[cal.obs[1]:(cal.obs[1]+diff(cal.obs)),]
    dates.m.c <- dates.gcm[cal.gcm[1]:(cal.gcm[1]+diff(cal.gcm)),]
    na.gcm <- rep(NA, n.gcm)

    # A very strange way of dividing up into equal chunks
    i.indices <- seq_along(lon)
    i.chunks <- suppressWarnings(matrix(c(i.indices, rep(NA, length(i.indices))),
                                        ncol=n.chunks*2))
    i.chunks <- i.chunks[,1:(min(which(is.na(i.chunks[1,])))-1), drop=FALSE]
    n.chunks <- ncol(i.chunks)

    jj <- seq_along(lat)

    browser()

    for(chunk in 1:2) { ##ncol(i.chunks)){
        ii <- na.omit(i.chunks[,chunk])
        ## I'm pretty sure that this chunking is backwards from optimal. The larger the ii the better
        ## Experiment measure 1068 x 10 vs 10 x 1068 (each are just under a GB)
        ## Actualy, they seem to be about the same, so, whatever
        cat('--> bias correcting', varname, 'chunk', chunk, '/', n.chunks, '-')
        o.c.chunk <- ncvar_get(obs, start=c(ii[1], jj[1], cal.obs[1]),
                               count=c(length(ii), length(jj), n.obs),
                               varid=varname, collapse_degen=FALSE)
        cat('-')
        m.p.chunk <- ncvar_get(gcm, start=c(ii[1], jj[1], 1),
                               count=c(length(ii), length(jj), n.gcm),
                               varid=varname, collapse_degen=FALSE)

        # FIXME: dont bother permuting and redimming... just do an apply(preserve=3) below
        m.p.chunk <- aperm(m.p.chunk, c(3, 1, 2))
        dim(m.p.chunk) <- c(n.gcm, length(ii) * length(jj))
        o.c.chunk <- aperm(o.c.chunk, c(3, 1, 2))
        dim(o.c.chunk) <- c(n.obs, length(ii) * length(jj))
        cat('*\n')
        ij <- t(expand.grid(i=seq_along(ii), j=seq_along(jj)))

        m.p.chunk <- foreach(ij=ij,
                             o.c=apply(o.c.chunk, 2, identity),
                             m.p=apply(m.p.chunk, 2, identity),
                             .inorder=TRUE,
                             .combine=cbind,
                             .multicombine=TRUE
                             #.options.mpi=mpi.options) %do% {
                             ) %do {
                                 i = ij['i',]
                                 j = ij['j',]
                                 if(all(is.na(o.c)) || all(is.na(m.p))) {
                                     na.gcm
                                 } else {
                                     m.c <- m.p[cal.gcm[1]:(cal.gcm[1]+diff(cal.gcm))]
                                     tQPQM(o.c=o.c, m.c=m.c, m.p=m.p, dates.o.c=dates.o.c,
                                           dates.m.c=dates.m.c, dates.m.p=dates.gcm,
                                           n.window=n.window, ratio=ratio[[varname]], trace=trace,
                                           jitter.factor=jitter.factor, seasonal=seasonal[[varname]],
                                           multiyear=multiyear, n.multiyear=n.multiyear,
                                           expand.multiyear=expand.multiyear, n.tau=pr.n.tau)
                                 }
                             }
        dim(m.p.chunk) <- c(n.gcm, length(ii), length(jj))
        cat(dim(m.p.chunk))
        m.p.chunk <- aperm(m.p.chunk, c(2, 3, 1))
        cat('--> writing chunk', chunk, '/', n.chunks, '\n')
        ncvar_put(nc=out, varid='pr', vals=m.p.chunk,
                  start=c(ii[1], jj[1], 1), count=dim(m.p.chunk))
        print(object.size(x=lapply(ls(), get)), units="Mb")
        gc()
    }

    nc_close(gcm)
    nc_close(obs)
    nc_close(out)

    print('Elapsed Time')
    print(proc.time() - ptm)
}
################################################################################
