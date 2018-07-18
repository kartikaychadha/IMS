# Remove_Overlaps.R
#
# Version:   0.1
# Author:    Kartikay Chadha
#
# Git rel-Date:    2018-07-17
# ====================================================================
#
#-------------------------CODE 1 (THIS IS NOT A FUNCTION)-----------------------------
#
# Call- name:   Remove_Overlaps(data,x)
# Purpose:      Remove SNPs with possible overlaps in flanking sequence. 
# Note:         Present in working directory
#               
#
# Depends:      randomRows()
#
# Inputs:       Dataframe (X, tag_rs, tag_chr, tag_pos, proxy_rs, proxy_chr, proxy_pos)
#               Order of columns CANNOT be altered.
#               x- flanking length (=BP 1 side x 2)  e.g: 201 (refer Thesis Chapter 3)        
# 
#
# Output:       Dataframe (list of all SNPs atleast x bp apart from each other)
#               Additional processing comlums (Flag values etc)
#       
#


Remove_overlaps <- function(data,x){
  
  f_list <- data.frame(X= numeric(),tag_rs =character(), tag_chr= numeric(), tag_pos= numeric(),proxy_rs =character(), proxy_chr= numeric(), proxy_pos= numeric(), dis = numeric(), flag = numeric() )
  m=1
  #Find unique tag SNps 
  unique_tag <- unique(data$tag_rs)
  
  #For each tag snp
  for(i in 1:length(unique_tag))
  {
    
    region <- subset(data,data$tag_rs==unique_tag[i])    #subset of data for each tag snp
    region <- region[ order(region$proxy_pos), ]         # Ordering in accesing order of proxy SNP
    region$dis[1] <- x                                   # Setting the starting SNP
    region$flag[1] <- 0
    ## termination for only one value 
    if(nrow(region)<2)
    {
      f_list<- rbind(f_list,region)
      
    }
    else
    {
      # Calculating intermarker distance
      for (j in 2:nrow(region))
      {
        region$dis[j] <- as.numeric(region$proxy_pos[j])-as.numeric(region$proxy_pos[j-1])
      }# intermarker distance calculated
      
      
      #Find clusters
      m=1
      region$flag=0
      for (j in 1:(nrow(region)-1))
      {
        #For all SNPs in the region
        if(region$dis[j] > x+1 && region$dis[j+1] > x+1)
        {
          region$flag[j]= 0
        }
        if(region$dis[j] < x+1 && region$dis[j+1] < x+1)
        {
          region$flag[j]= m
        }
        if(region$dis[j] > x+1 && region$dis[j+1] < x+1)
        {
          region$flag[j]= m
        }
        if(region$dis[j] < x+1 && region$dis[j+1] > x+1)
        {
          region$flag[j]= m
          m= m+1
        }
        
        #for the last SNP
        if(j==nrow(region)-1)
        {
          if(region$dis[j+1] < x+1)
          {
            region$flag[j+1]= m
          }
          else
          {
            region$flag[j+1]= 0
          }
        }
      }#for cluster loop ends
      
      
      
      #update f-list and region list
      f_list<- rbind(f_list,subset(region,region$flag==0))
      region<-subset(region,region$flag!=0)
      
      #for each cluster in a region 
      unique_clusters <- unique(region$flag)
      for(j in 1:length(unique_clusters))
      {
        subcluster <- data.frame()
        subcluster <- subset(region,region$flag == j)
        while(nrow(subcluster)!=0)
        {
          select <- as.character(sample(subcluster$proxy_rs,1))
          f_list<- rbind(f_list,subcluster[subcluster$proxy_rs==select,])
          
          #for lower values
          k=1
          while(k<=nrow(subcluster))
          {
            if(subcluster$proxy_rs[k] == select || is.na(subcluster$proxy_rs[k])){break}
            if(subcluster$proxy_pos[subcluster$proxy_rs==select]-subcluster$proxy_pos[k] < x+1)
            {
              subcluster <- subset(subcluster,subcluster$proxy_rs!=subcluster$proxy_rs[k])
              k=k-1
            }
            k=k+1
          }
          #for upper values
          k=nrow(subcluster)
          while(k>=1)
          {
            if(subcluster$proxy_rs[k] == select || is.na(subcluster$proxy_rs[k])){break}
            if(subcluster$proxy_pos[k]-subcluster$proxy_pos[subcluster$proxy_rs==select] < x+1)
            {
              subcluster <- subset(subcluster,subcluster$proxy_rs!=subcluster$proxy_rs[k])
              k=k+1
            }
            k=k-1
          }
          subcluster <- subset(subcluster,subcluster$proxy_rs!=select)
        }
        
        
      }
    } 
    message(i,"\r\n",appendLF=FALSE)
    flush.console()
  }
  
  
  return (f_list)
  
}