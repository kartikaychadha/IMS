# Extract sequences from csv file- for GTEX permuations
#
#

#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

#Change number of runs here 
foreach(i=1:1000,.export = 'fread') %dopar% {
  write.csv(fread(paste(1," _AOI.csv",sep = ""),select = c(6)),paste(i,"SOI.csv",sep = ""),row.names = F)
  write.csv(fread(paste(1," _ANOI.csv",sep = ""),select = c(6)),paste(i,"CS.csv",sep = ""),row.names = F)
}
  

#stop cluster
stopCluster(cl)