# Extract sequences from csv file- for GTEX permuations
#
#

#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

#Change number of runs here 
#foreach(i=1:1000,.export = 'fread') %dopar% {
for(i in 1:1000){
  write.csv(fread(paste(i," _AOI.csv",sep = ""),select = c(6)),paste(i,"SOI.csv",sep = ""),row.names = F)
  write.csv(fread(paste(i," _ANOI.csv",sep = ""),select = c(6)),paste(i,"CS.csv",sep = ""),row.names = F)
print(i)
  }
  

#stop cluster
stopCluster(cl)