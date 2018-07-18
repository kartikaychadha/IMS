# dependencies.R
#
# Purpose:   Incudes all functions may be needed to perform DNA word analysis. 
# Version:   0.1
# Author:    Kartikay Chadha
#
# Input:     Function specific- Read below 
# Output:    Function specific- Read below 
# Depends:   N/A
#
# ToDo:      Run all functions in this file before performing running the analysis. 
# Notes:     
#
# Git rel-Date:    2018-07-17
# ====================================================================
#
#-------------------------Function 1-----------------------------
#
# Call- name:   find_overlaps(pattern,string)
# Purpose:      To claculate overlaps
#
# Desc:         Search occurance of p in s; including overlaps 
#               Eg: Search p= "AAAA" in s="AAAAAAAT" O/P: 4
#
# Depends:      NA
# Inputs:       p (pattern) to be searched in s (string)- p, s 
# Output:       Count (numeric)


find_overlaps <- function(p,s) {
  gg <- gregexpr(paste0("(?=",p,")"),s,perl=TRUE)[[1]]
  if (length(gg)==1 && gg==-1) 0 else length(gg)
}

#
#
#-------------------------Function 2-----------------------------
#
# Call- name:   CochranArmitageTest(Contigency Table, paramters)
# Purpose:      Cochran Armitage Trend Test
#
# Desc:         Manual Link:
#               http://cran.wustl.edu/web/packages/DescTools/DescTools.pdf
#
# Depends:      NA
# Inputs:       Numeric Table (Contigency table)
#               parameters (two.sided,increasing,decreasing)
# Output:       Statictics, DIM, p-value
#
#Citation:      Andri Signorell et mult. al. (2016). DescTools: Tools for descriptive statistics. R package version 0.99.17
#

CochranArmitageTest <- function(x, alternative = c("two.sided","increasing","decreasing")) {
  
  # based on:
  # http://tolstoy.newcastle.edu.au/R/help/05/07/9442.html
  DNAME <- deparse(substitute(x))
  
  if (!(any(dim(x)==2)))
    stop("Cochran-Armitage test for trend must be used with rx2-table", call.=FALSE)
  
  if (dim(x)[2]!=2) x <- t(x)
  
  nidot <- apply(x, 1, sum)
  n <- sum(nidot)
  
  # Ri <- scores(x, 1, "table")
  Ri <- 1:dim(x)[1]
  Rbar <- sum(nidot*Ri)/n
  
  s2 <- sum(nidot*(Ri-Rbar)^2)
  pdot1 <- sum(x[,1])/n
  z <- sum(x[,1]*(Ri-Rbar))/sqrt(pdot1*(1-pdot1)*s2)
  STATISTIC <- z
  
  alternative <- match.arg(alternative)
  
  PVAL <- switch(alternative,
                 two.sided = 2*pnorm(abs(z), lower.tail=FALSE),
                 increasing = pnorm(z),
                 decreasing = pnorm(z, lower.tail=FALSE) )
  
  PARAMETER <- dim(x)[1]
  names(STATISTIC) <- "Z"
  names(PARAMETER) <- "dim"
  
  METHOD <- "Cochran-Armitage test for trend"
  structure(list(statistic = STATISTIC, parameter = PARAMETER, alternative = alternative,
                 p.value = PVAL, method = METHOD, data.name = DNAME
  ), class = "htest")
  
}

#
#
#-------------------------Function 3-----------------------------
#
# Call- name:   randomRows(DataFrame, NumberOfsequences)
# Purpose:      Randomly select rows in a dataFrame
#
# Depends:      NA
#
# Inputs:       Dataframe with multiple rows
#               Number of Rows to be returned 
#
# Output:       Random rows returned 
#


randomRows<-function(df,n){
  return(df[sample(nrow(df),n),])
}
