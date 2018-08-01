# FishersOnData.R
#
# Version:   0.1
# Author:    Kartikay Chadha
#
# Git rel-Date:    2018-08-01
# ====================================================================
#
#-------------------------Function 1-----------------------------
#
# Call- name:   FishersOnData(Sequences of Interest (SOI) count Matrix, Control Seq count matrix, Total groups of DNA words)
# Purpose:      Pereforms ONLY Fisher's Exact test.
#
# Note:         
#
# Depends:      Proportion_DNAwords.R
# Inputs:       Output from Count_DNAwords.R
#               1. Matrix for SOI from Count_DNAwords.R
#               2. Matrix for Control Sequences from Count_DNAwords.R
#
#               Dataframe structure: (each cell shows the frequency corrusponding to group of DNA words and DNA sequence)
#               Columns : Groups of DNAwords
#               Rows    : DNA seqeunces 
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
FishersOnData <- function(matrix_AOI, matrix_ANOI, no_of_features ) {
  
  chi_fisher_cells <- data.frame(feature_no = numeric(0),VOI_pre = numeric(0),VOI_abs = numeric(0), CS_pre =numeric(0),CS_abs = numeric(0), chi_fis = character(0), stat_OR = numeric(0), p_value = numeric(0), stringsAsFactors = F)
  for(i in 1:no_of_features)
  {
    chi_fisher_cells[i,1]<-i
    f1<-Proportion_DNAwords(matrix_AOI,i)
    f2<-ncol(matrix_AOI)-f1
    f3<-Proportion_DNAwords(matrix_ANOI,i)
    f4<-ncol(matrix_ANOI)-f3
    chi_fisher_cells[i,2]<-f1
    chi_fisher_cells[i,3]<-f2
    chi_fisher_cells[i,4]<-f3
    chi_fisher_cells[i,5]<-f4
    
    chi_fisher_cells[i,6]<-"fis"
      
      c=fisher.test(rbind(c(f1,f2),c(f3,f4)), alternative="two.sided",workspace=1e19)$estimate
      chi_fisher_cells[i,7]<-c
       remove(c)
       d= fisher.test(rbind(c(f1,f2),c(f3,f4)), alternative="two.sided",workspace=1e19)$p.value
      
      
      chi_fisher_cells[i,8]<-d
      
    remove(d)
    
    
  }
  return(chi_fisher_cells)
}