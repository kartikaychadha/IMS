#Freq_DNAwords.R
#
# Version:   0.1
# Author:    Kartikay Chadha
#
# Git rel-Date:    2018-07-17
# ====================================================================
#
#-------------------------Function 1-----------------------------
#
# Call- name:   Freq_DNAwords(fasta_data_frame, DNAword, Threshold freq cut-off)
# Purpose:      Returns the frequency of a DNA word in all Fasta sequences. 
#
#
# ToDo:         Define minimum frequecy as m
# Note:         minimum frequency = 1 (default) i.e., atleast 1 time occurance in a sequence. 
#               To run code for all sequences: define m = -1
#
# Depends:      find_overlaps() --Function 1 
#
# Inputs:       data - Fasta file (data frame) [Columns=1; Rows=Number of seq x2]-- SEE Example Files/fasta_example.txt
#               DNAword : Single string entry. Eg: "AAAAA" 
#               m : Threshold cut-off frequency 
#
#
# Outputs:      Dataframe
#               Column 1: description line of the seqeunce
#               Column 2: Frequency of DNA words in that sequence 
#
#

Freq_DNAwords <- function(data, DNAword, m = 1){
  k <- 0
  x <- data.frame()
  counter <- 1
  for(i in 1:(nrow(data)/2))
  {
    k=k+2
    if(find_overlaps(as.character(DNAword),data[k,1])>m)
    {
      x[counter,1]<-data[k-1,1]
      x[counter,2]<-find_overlaps(as.character(DNAword),data[k,1])
      counter = counter + 1 
    }
    
  }
  return (x)
}
