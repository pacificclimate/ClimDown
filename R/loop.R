cores <- getOption('mc.cores')
if (cores > 1) {
  `%loop%` <- `%dopar%`
  start.par.backend <- function() {
    backend <-  getOption('par.backend')
    if ( backend == 'MPI') {
      require(doMPI)
      cl <- startMPIcluster(cores-1)
      registerDoMPI(cl)
    }
    else if (backend == 'multicore') {
      require(doParallel)
      registerDoParallel(cores=cores)
    } else {
      stop('Option "par.backend" must be "multicore" or "MPI"')
    }
    cl
  }
  stop.par.backend <- function(cl) {
    backend <-  getOption('par.backend')
    if ( backend == 'MPI') {
      require(doMPI)
      closeCluster(cl)
    } else if (backend == 'multicore') {
      require(doParallel)
      stopImplicitCluster()
    } else {
      stop('Option "par.backend" must be "multicore" or "MPI"')
    }
  }
} else {
  `%loop%` <- `%do%`
  start.par.backend <- function (){}
  stop.par.backend <- function(cl) {}
}
