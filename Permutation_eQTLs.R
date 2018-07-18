# Permutation_eQTLs.R
#
# Version:   0.1
# Author:    Kartikay Chadha
#
# Git rel-Date:    2018-07-17
# ====================================================================
#
#-----------------------CODE 1 (THIS IS NOT A FUNCTION)------------------
#
# Call- name:   Remove_Overlaps(data,x)
# Purpose:      Remove SNPs with possible overlaps in flanking sequence. 
# Note:         This code is designed for 9,103 Thy eQTL. Please read and change code as per requirement. 
#               Present in working directory.
#
# Depends:      randomRows()
#
# Inputs:       Dataframe for ALL Variants of Interest (VOI) and Control varinats. 
#               Order of columns CANNOT be altered.
#               Colums: Chr_position, rsID, Allele1, MAF, Allele2, DNA sequence 
#                       
# 
#
# Output:       Dataframes: Files from each run (eg: 1_AOI.csv and 2_ANOI.csv and 1_res.csv)
#              
#       
#


#Load the VOI and control list together
AOI_ANOI_together_sorted<-fread("../AOI_ANOI_together_sorted.txt")
#Load the list of groupd of DNA words
features<-fread("../features.txt")

#1-1000 runs
#CHANGE NUMBER OF RUNS HERE
foreach(i=1:1000) %dopar% {
  
  # Randomly shuffle all VOI and control variants
  AOI_ANOI_together_sorted<-AOI_ANOI_together_sorted[sample(1:nrow(AOI_ANOI_together_sorted)), ]
  
  #Split the shuffled file into two equal half of 9103
  #Change "9,103" HERE AS PER YOUR REQUIREMENTS
  
  x1<-AOI_ANOI_together_sorted[1:9103,]
  x2<-AOI_ANOI_together_sorted[9104:18206,]
  
  #Documenting the new VOI and control variants
  write.csv(x1,paste(i,"_AOI.csv"),row.names = F)
  write.csv(x2,paste(i,"_ANOI.csv"),row.names = F)
  
  res<-data.frame(stringsAsFactors = F)
  
  foreach(j=1:nrow(features)) %dopar% {
    
    #SOI counting   
    m1_AOI<-length(grep(features[j,1],x1$V6))
    
    if(!is.na(features[j,2])){
      m2_AOI<-length(grep(features[j,2],x1$V6))
    }else
    {
      m2_AOI=0
    }
    
    #Control Seq counting   
    m1_ANOI<-length(grep(features[j,1],x2$V6))
    
    if(!is.na(features[j,2])){
      m2_ANOI<-length(grep(features[j,2],x2$V6))
    }else
    {
      m2_ANOI=0
    }
    
    
    
    
    f1<-m1_AOI+m2_AOI
    f2<-9103-f1
    f3<-m1_ANOI+m2_ANOI
    f4<-9103-f3
    
    
    res[j,1]<-features[j,1]
    res[j,2]<-features[j,2]
    res[j,3]<-f1
    res[j,4]<-f2
    res[j,5]<-f3
    res[j,6]<-f4
    res[j,7]<-"fis"
    res[j,8]<-fisher.test(rbind(c(f1,f2),c(f3,f4)))$estimate
    res[j,9]<-fisher.test(rbind(c(f1,f2),c(f3,f4)))$p.value
    
    print(j)
    
  }
  write.csv(res,paste(i,"_res.csv"),row.names = F)
  remove(rm)
  remove(x1)
  remove(x2)
  print(paste("**********************",i,sep = "      "))
} 