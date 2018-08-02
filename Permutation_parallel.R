# Permutation_parallel.R
#
# Version:   0.1
# Author:    Kartikay Chadha
#
# Git rel-Date:    2018-08-02
# ====================================================================
#
#-------------------------THIS IS NOT A FUNCTION-----------------------------
#
# Call- name:   NA
# Purpose:      Counts the words in the sequences. Performs fisher exact test and reports OR + P-value. 
#
# Note:               
#
# Depends:      FisherUsingFasta()
# 
# To-do:        Set working directory with all fasta files as .csv
#               
# Inputs:       SOI - Fasta files (.csv) [Columns=1; Rows=Number of seq x2]-- SEE Example Files/fasta_example.txt 
#               Control sequnces - Fasta file (.csv) [Columns=1; Rows=Number of seq x2]-- SEE Example Files/fasta_example.txt 
#               groups of DNA words (data frame)[Columns=2; Rows=Number of possible sequences]--See Example Files/groups_DNAwords_5mer.csv
#               n=number of permuations
#              
#
#
# Output:       Dataframes X n  
#               Columns : DNAwordGroupID, 
#                         #SOI containing groupDNAword 
#                         #SOI NOT containing groupDNAword
#                         #Control sequences containing groupDNAword
#                         #Control sequences NOT containing groupDNAword
#                         Test Statistic
#                         P-value 
#
#               Rows    : Groups of DNA words 
#




library(foreach)
library(doParallel)


#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

#Change number of runs here 
foreach(i=1:1000,.export = 'fread') %dopar% {
  
  #Reading files 
  #Change file names here
  m<-fread(paste(i,"AOI.csv",sep = ""),select = c(2))
  n<-fread(paste(i,"ANOI.csv", sep = ""),select = c(2))
  
  #Performing analysis 
  result<-FisherUsingFasta(m,n,groups_DNAwords_5mer)
  write.csv(result,paste(i,"Result.csv",sep = ""),row.names = F)
  remove(m)
  remove(n)
  
 
  
}

#stop cluster
stopCluster(cl)
