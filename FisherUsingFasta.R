# FisherUsingFasta.R
#
# Version:   0.1
# Author:    Kartikay Chadha
#
# Git rel-Date:    2018-08-01
# ====================================================================
#
#-------------------------Function 1-----------------------------
#
# Call- name:   Analysis_short(Sequences of Interest (SOI) fasta, Control Seq count fasta, Groups of DNA words)
# Purpose:      Counts the words in the sequences. Performs fisher exact test and reports OR + P-value. 
#
# Note:               
#
# Depends:      NA
#               
# Inputs:       SOI - Fasta file (data frame) [Columns=1; Rows=Number of seq x2]-- SEE Example Files/fasta_example.txt 
#               Control sequnces - Fasta file (data frame) [Columns=1; Rows=Number of seq x2]-- SEE Example Files/fasta_example.txt 
#               groups of DNA words (data frame)[Columns=2; Rows=Number of possible sequences]--See Example Files/groups_DNAwords_5mer.csv
# 
#              
#
#
# Output:       Dataframe
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


FisherUsingFasta <- function(SOI_fasta, CS_fasta, Groups){
  
  #Removing descrption lines 
  SOI_fasta<-SOI_fasta[seq(2, nrow(SOI_fasta), 2),]
  CS_fasta<-CS_fasta[seq(2, nrow(CS_fasta), 2),]
  
  fisher <- data.frame(feature_no = numeric(0),VOI_pre = numeric(0),VOI_abs = numeric(0), CS_pre =numeric(0),CS_abs = numeric(0),OR = numeric(0), p_value = numeric(0), stringsAsFactors = F)
  
  for(j in 1:nrow(Groups)){
    fisher[j,1]<-j
    fisher[j,2]<-sum(mapply(grepl,paste(Groups[j,1],Groups[j,2],sep = "|"),SOI_fasta), na.rm=TRUE)
    fisher[j,3]<-nrow(SOI_fasta)-fisher[j,2]
    fisher[j,4]<-sum(mapply(grepl,paste(Groups[j,1],Groups[j,2],sep = "|"),CS_fasta), na.rm=TRUE)
    fisher[j,5]<-nrow(CS_fasta)-fisher[j,4]
    
    #fisher[j,6]<-fisher.test(rbind(c(fisher[j,2],fisher[j,3]),c(fisher[j,4],fisher[j,5])), alternative="two.sided",workspace=1e19)$estimate
   #fisher[j,7]<-fisher.test(rbind(c(fisher[j,2],fisher[j,3]),c(fisher[j,4],fisher[j,5])), alternative="two.sided",workspace=1e19)$p.value
    
    
  }
  fisher$OR<-apply(fisher,1,function(x) fisher.test(matrix(x[2:5],nrow=2))$estimate)
  fisher$p_value<-apply(fisher,1,function(x) fisher.test(matrix(x[2:5],nrow=2))$p.value)  
  fisher<-cbind(Groups,fisher)
  return(fisher)
  
}
