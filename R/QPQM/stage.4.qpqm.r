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

    nlat <- qpqm$dim$lat$len
    nlon <- qpqm$dim$lon$len

    bcca.time <- compute.time.stats(nc)

    n.qpqm <- qpqm$dim$time$len
    n.bcca <- bcca.time$n

    years <- unique(format(bcca.time$vals, '%Y'))

    for (year in years) { # One year at a time
        print(year)
        i <- which(year == format(bcca.time$vals, '%Y')
        t.st <- i[1]
        print(t.st)
        date.sub <- bcca.time$vals[i]
        t.cnt <- length(date.sub)
        mons <- as.numeric(format(date.sub, '%m'))

        read.time <- proc.time()
        # Read one year for all space
        var.bcca <- ncvar_get(bcca,varid=var.name,start=c(1,1,t.st),count=c(nlon,nlat,t.cnt))
        var.qpqm <- ncvar_get(qpqm,varid=var.name,start=c(1,1,t.st),count=c(nlon,nlat,t.cnt))
        var.final <- var.qpqm*0

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

        par.time <- proc.time()

        dqm <- foreach(i=1:12,
                       b.list=bcca.list,
                       q.list=qpqm.list,
                       .inorder=TRUE,
                       .options.mpi=mpi.options) %dopar% {

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
