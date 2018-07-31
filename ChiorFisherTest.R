# ChiorFisherTest.R
#
# Version:   0.1
# Author:    Kartikay Chadha
#
# Git rel-Date:    2018-07-17
# ====================================================================
#
#-------------------------Function 1-----------------------------
#
# Call- name:   ChiorFisherTest(Sequences of Interest (SOI) count Matrix, Control Seq count matrix, Total groups of DNA words)
# Purpose:      Pereforms Fisher's Exact test or Chi-Squared test
#
# Note:         Any cell's value less than 5 triggers Fishers exact test else Chi squared is performed in this function. 
#               
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
ChiorFisherTest <- function(matrix_AOI, matrix_ANOI, no_of_features ) {
  
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
    
    if(f1<5 || f2<5 || f3<5 || f4<5)
    {
      
      c=fisher.test(rbind(c(f1,f2),c(f3,f4)), alternative="two.sided")$estimate
      d= fisher.test(rbind(c(f1,f2),c(f3,f4)), alternative="two.sided")$p.value
      chi_fisher_cells[i,6]<-"fis"
      chi_fisher_cells[i,7]<-c
      chi_fisher_cells[i,8]<-d
      
    }
    else{
      c=chisq.test(rbind(c(f1,f2),c(f3,f4)), simulate.p.value=F)$statistic
      d=chisq.test(rbind(c(f1,f2),c(f3,f4)), simulate.p.value=F)$p.value
      chi_fisher_cells[i,6]<-"chi"
      chi_fisher_cells[i,7]<-c
      chi_fisher_cells[i,8]<-d
    }
    
    
  }
  return(chi_fisher_cells)
}