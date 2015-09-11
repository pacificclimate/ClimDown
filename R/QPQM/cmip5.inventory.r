###Script to assemble a human readable inventory of CMIP5 data at PCIC

base.dir <- '/home/data/climate/CMIP5/'

centres <- list.files(path=base.dir)

if (1==1) {
full.list <- c('Centre','Model','Scenario','Run','Variable','Start','End')
variable <- 'tasmin'
scen.list <- c('historical','rcp26','rcp45','rcp60','rcp85')


for (c in seq_along(centres)) {
  centre <- centres[c]
  print(centre)
  ##PR files
  pr.files <- as.list(list.files(path=paste(base.dir,centre,sep=''),pattern=paste(variable,'_day',sep=''),recursive=TRUE))
  if (length(pr.files)!=0) {
    print(pr.files)
    fsplit <- lapply(pr.files,function(x) {strsplit(x,'/')[[1]]})
    file.matrix <- matrix(NA,nrow=length(pr.files),ncol=length(fsplit[[1]])-1)
    for (i in 1:length(pr.files))
      file.matrix[i,] <- fsplit[[i]][1:(length(fsplit[[1]])-1)]
    
    nc.files <- unlist(lapply(pr.files,function(x) {strsplit(x,'/')[[1]][9]}))
    
    years <- unlist(regmatches(nc.files,gregexpr('[0-9]{8}-[0-9]{8}',nc.files)))
    yst <- substr(years,1,4)
    yen <- substr(years,10,13)
    
    file.matrix <- cbind(rep(centre,length(pr.files)),file.matrix,yst,yen)
    file.sub <- file.matrix[,c(1,2,3,7,9,10,11)]
    ##names(file.sub) <- c('Centre','Model','Scenario','Run','Variable','Start','End')
    
    scen.sub <- file.sub[,3] %in% scen.list
    
    rv.matrix <- file.sub[scen.sub,]
    full.list <- rbind(full.list,rv.matrix)
  }
}
}
##Sort the data into a more readable format
read.list <- c('Model','Scenario','Runs','Years')
models <- unique(full.list[-1,2])

for (m in seq_along(models)) {
  model <- models[m]
  sub.list <- full.list[full.list[,2]==model,]
  if (is.null(dim(sub.list))) {
  } else {
    scenarios <- unique(sub.list[,3])
    for (scen in scenarios) {
      scen.list <- sub.list[sub.list[,3]==scen,]
      if (is.null(dim(scen.list))) {
        runs <- scen.list[4]
        years <- paste(scen.list[6],scen.list[7],sep='-')
      } else {
        runs <- paste(scen.list[,4],collapse='|')
        bnds <- cbind(as.character(scen.list[,6]),as.character(scen.list[,7]))
        years <- apply(bnds,1,paste,collapse='-')
      }
      if (length(unique(years))>1)
        warning(paste('Years in models runs are not the same for: ', model,'-',scen,': ',years,sep=''))
      rv <- as.character(c(model,scen,runs,years[1]))
      read.list <- rbind(read.list,rv)
    }
  }
}

write.table(read.list,file=paste(variable,'_CMIP5_model_inventory.csv',sep=''),quote=FALSE,row.name=F,col.name=F,sep=',')
