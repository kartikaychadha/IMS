# Proportion_DNAwords.R
#
# Version:   0.1
# Author:    Kartikay Chadha
#
# Git rel-Date:    2018-07-17
# ====================================================================
#
#-------------------------Function 1-----------------------------
#
# Call- name:   Proportion_DNAwords(data, DNA_group_ID)
# Purpose:      Returns the proportion of sequnces matching a specific groups of DNA Word. 
#
#
# ToDo:         Run this code for all groups of DNAwords i.e., 1 to ncol(data)      
# Note:         DNAword ID can be extracted fromt DNA:
#               groups of DNA words (data frame)[Columns=2; Rows=Number of possible sequences]--See Example Files/groups_DNAwords_5mer.csv
#
# Depends:      Data output from Count_DNAwords.R
#
# Inputs:       Dataframe (OUTPUT FROM Count_DNAwords.R)
#               Columns : Groups of DNAwords
#               Rows    : DNA seqeunces 
#
#               DNA_group_ID: Column number for the group of DNA word 
#
#
# Outputs:      Number of Sequences coninting the selected group of DNA words. 
#
#


Proportion_DNAwords <- function(data,feature_num){
  sum<-0
  for(i in 1:ncol(data))
  {
    if(data[feature_num,i]>0)
    {
      sum=sum+1;  
    }
  }
  
  return (sum)
}