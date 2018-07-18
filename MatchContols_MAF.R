# MatchControls_MAF.R
#
# Version:   0.1
# Author:    Kartikay Chadha
#
# Git rel-Date:    2018-07-17
# ====================================================================
#
#-------------------------CODE 1 (THIS IS NOT A FUNCTION)-----------------------------
#
# Call- name:   MatchControls_MAF
# Purpose:      Matches the Variants of interest (VOI) with control variants for only Minor Allele Frequency. 
# Note:         Present in working directory
#               
#
# Depends:      randomRows()
# Inputs:       Dataframe- List of VOI (Columns: Chr_position, rsID, Allele1, MAF, Allele2)
#               Dataframe- Pool of controls within which the Matching is performed (Columns: Chr_position, rsID, Allele1, MAF, Allele2)           
#               Order of column CANNOT be altered. 
# 
#
# Output:       Dataframe
#               Chr_position_VOI, rsID_VOI, Allele1_VOI, MAF_VOI, Allele2_VOI, Chr_position_Controls, rsID_Controls, Allele1_Controls, MAF_Controls, Allele2_Controls
#




# Reading files 
l <- fread("eqtl_final.txt")
k <- fread("Thy_p0.5_all_gtex_maf.txt")

# Creating result df
df <- data.frame(stringsAsFactors = F)

# for all VOI


for(i in 1:nrow(l))
{
  #Select the top most row from the VOI dataframe
  sam<-l[1,]
  while (length(grep(sam[1,2],l$V2))!=0){
    #Select a random row matching the criteria for MAF (0.01 +/-)
    #Change Widow in the next line.
    sam <- randomRows(subset(k,k$V4<(as.numeric(l[i,4])+0.01) & k$V4>(as.numeric(l[i,4])-0.01)),1)}
  df[i,1] <- l[i,1]
  df[i,2] <- l[i,2]
  df[i,3] <- l[i,3]
  df[i,4] <- l[i,4]
  df[i,5] <- l[i,5]
  
  df[i,6] <- sam[1,1]
  df[i,7] <- sam[1,2]
  df[i,8] <- sam[1,3]
  df[i,9] <- sam[1,4]
  df[i,10] <-sam[1,5]
  rm(sam)
  message(i,"\r\n",appendLF=FALSE)
  flush.console()
}
