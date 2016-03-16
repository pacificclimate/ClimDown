library(ncdf4)
library(abind)
library(PCICt)

reorder <- function(x,ix) {
  rx <- x[ix]
  return(rx)
}

# Quantile perturbation quantile mapping
qdm.netcdf.wrapper <- function(qpqm.file, bcca.file, analogues, out.file, varname='tasmax') {
    ## Arguments:
    ## qpqm.file - Gridded historical observations (i.e. QPQM output)
    ## gcm_file - GCM simulations interpolated to the obs_file grid (i.e. BCCA output)
    ## 
    ## output_file - The file to create (or overwrite) with the bias corrected outputs
    ## var.name - One of pr, tasmax, tasmin
    ##
    ## All files should have the same spatial domain.\n\n')
    ptm <- proc.time()

    qpqm <- nc_open(qpqm.file)
    bcca <- nc_open(bcca.file)
    out  <- nc_open(out.file,write=TRUE)

    lat <- qpqm$dim$lat$vals
    nlat <- length(lat)
    lon <- qpqm$dim$lon$vals
    nlon <- length(lon)

    dates.qpqm <- as.PCICt(strsplit(qpqm$dim$time$units, ' ')[[1]][3],
                           qpqm$dim$time$calendar) + qpqm$dim$time$vals*86400
    dates.qpqm <-  apply(do.call(rbind, strsplit(format(dates.qpqm, '%Y %m %d'),
                                                 ' ')), 2, as.integer)

    dates.bcca <- as.PCICt(strsplit(bcca$dim$time$units, ' ')[[1]][3],
                           bcca$dim$time$calendar) + bcca$dim$time$vals*86400
    dates.bcca <-  apply(do.call(rbind, strsplit(format(dates.bcca, '%Y %m %d'),
                                                 ' ')), 2, as.integer)

    n.qpqm <- nrow(dates.qpqm)
    n.bcca <- nrow(dates.bcca)

    years <- unique(dates.bcca[,1])
    n.chunks <- 2 ##length(years)

    ##----------------------------------------------------------
    nlon <- 400

    for (n in 1:n.chunks) { # One year at a time
        year <- years[n]
        print(year)
        t.st <- grep(year,dates.bcca[,1])[1]
        print(t.st)
        date.sub <- dates.bcca[dates.bcca[,1] %in% year,]
        t.cnt <- nrow(date.sub)
        mons <- date.sub[,2]

        read.time <- proc.time()
        # Read one year for all space
        var.bcca <- ncvar_get(bcca,varid=var.name,start=c(1,1,t.st),count=c(nlon,nlat,t.cnt))
        var.qpqm <- ncvar_get(qpqm,varid=var.name,start=c(1,1,t.st),count=c(nlon,nlat,t.cnt))
        var.final <- var.qpqm*0

        # This does nothing and is a waste of time
        ncvar_put(nc=out, varid=var.name, vals=var.final,
                  start=c(1, 1, t.st), count=dim(var.final))
        print('Reading data time')
        print(proc.time()-read.time)

        print('Arrange into list')
        mn.subset <- vector(mode='list',length=12)
        bcca.list <- vector(mode='list',length=12)
        qpqm.list <- vector(mode='list',length=12)
        # Split the year up into 12 months
        for (m in 1:12) {
            subset <- which(mons==m)
            mn.subset[[m]] <- subset
            bcca.list[[m]] <- var.bcca[,,subset]
            qpqm.list[[m]] <- var.qpqm[,,subset]
        }
        ##qpqm.list=qpqm.list,
        print('Parallel')
        par.time <- proc.time()
        ##for (i in 1:12) {
        dqm <- foreach(i=1:12,
                       b.list=bcca.list,
                       q.list=qpqm.list, ##, ##) %dopar% {
                       .inorder=TRUE,
                       .options.mpi=mpi.options) %dopar% {
                           ##.combine='cbind',
                           ##.multicombine=TRUE,

                           b.list <- jitter(bcca.list[[i]],0.01) # Why the jitter on the BCCA output?
                           q.list <- qpqm.list[[i]]

                           slen <- dim(b.list)[3]
                           nlon <- dim(b.list)[1]
                           nlat <- dim(b.list)[2]
                           # FIXME: We need to actually apply the BCCA analogues to the obs here
                           # For each point in space, rank the days within a single month
                           bcca.ranks <- round(aperm(apply(b.list,c(1,2),rank,ties.method='average'),c(2,3,1))) ##Permute to fix the rank permute
                           data.p <- aperm(q.list,c(3,1,2))

                           # Permute the QPQM data into time x cell...
                           dim(data.p) <- c(slen,nlon*nlat)
                           # ... then convert it to a list
                           data.l <- lapply(apply(data.p,2,as.list),unlist)

                           # Permute the BCCA ranks as well into time x cell...
                           index.p <- aperm(bcca.ranks,c(3,1,2))
                           dim(index.p) <- c(slen,nlon*nlat)
                           # ... and make them a list
                           index.l <- lapply(apply(index.p,2,as.list),unlist)

                           # For each cell, reorder the time steps based on the ranks
                           result <- mapply(reorder,data.l,index.l)
                           dim(result) <- c(slen,nlon,nlat)
                           final <- aperm(result,c(2,3,1))
                       }
        for (l in 1:12)
            var.final[,,mn.subset[[l]]] <- dqm[[l]]

        ncvar_put(nc=out, varid=var.name, vals=var.final,
                  start=c(1, 1, t.st), count=dim(var.final))
    }

    nc_close(qpqm)
    nc_close(bcca)
    nc_close(out)

    print('Elapsed time')
    print(proc.time() - ptm)
}
