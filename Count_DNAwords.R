# Count_DNAwords.R
#
# Version:   0.1
# Author:    Kartikay Chadha
#
# Git rel-Date:    2018-07-17
# ====================================================================
#
#-------------------------Function 1-----------------------------
#
# Call- name:   Count_DNAwords(fatsa_data_frame, groups_of_DNAwords)
# Purpose:      Count the occurance of all groups of DNA words in all fasta seqeunces
#
# Note:         Count the occurance of each group of DNA word in each sequence and returns a matrix:
#               Columns : Groups of DNAwords
#               Rows    : DNA seqeunces 
#               
#
# Depends:      find_overlaps()
# Inputs:       data - Fasta file (data frame) [Columns=1; Rows=Number of seq x2]-- SEE Example Files/fasta_example.txt 
#               groups of DNA words (data frame)[Columns=2; Rows=Number of possible sequences]--See Example Files/groups_DNAwords_5mer.csv
#
#
# Output:       Dataframe (each cell shows the frequency corrusponding to group of DNA words and DNA sequence)
#               Columns : Groups of DNAwords
#               Rows    : DNA seqeunces 
#


Count_DNAwords <- function(fasta_data_frame, features_data_frame){
  
  k<-0
  result<-data.frame()
  for(i in 1:(nrow(fasta_data_frame)/2))
  {
    k<-k+2
    
    for(j in 1:nrow(features_data_frame))
    {
      count<-0
      count<-find_overlaps(features_data_frame[j,1],fasta_data_frame[k,1])+find_overlaps(features_data_frame[j,2],fasta_data_frame[k,1])
      result[j,i]<-count
      
    }
    
  }
  return(result)
} 