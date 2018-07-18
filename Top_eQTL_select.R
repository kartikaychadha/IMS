# Top_eQTL_select.R
#
# Version:   0.1
# Author:    Kartikay Chadha
#
# Git rel-Date:    2018-07-17
# ====================================================================
#
#-------------------------CODE 1 (THIS IS NOT A FUNCTION)-----------------------------
#
# Call- name:   Top_eQTL_select
# Purpose:      Returns the SNP-gene pair with minimum P-value. Only one pair is returned incase more than 1 pair has Pval = Min(Pval)
#
# Note:         
#               
#
# Depends:      randomRows()
# Inputs:       GTEX FILE as dataframe named: Thyroid_Analysis_v6p_signif_snpgene_pairs
#               (Eg: javascript:portalClient.browseDatasets.downloadFile('GTEx_Analysis_v6p_eQTL.tar','gtex_analysis_v6p/single_tissue_eqtl_data/GTEx_Analysis_v6p_eQTL.tar'))
#
#
# Output:       Dataframe (same format as input)
#               
#

#result df
Top_eqtls_thyroid<-data.frame(stringsAsFactors = F)

k=1
#code
for (i in unique(Thyroid_Analysis_v6p_signif_snpgene_pairs$gene_id))
{
  genewise <- subset(Thyroid_Analysis_v6p_signif_snpgene_pairs,Thyroid_Analysis_v6p_signif_snpgene_pairs$gene_id==i)
  genewise_min <- subset(genewise,genewise$pval_nominal==min(genewise$pval_nominal))
  
  if(nrow(genewise_min)>1)
  {
    Top_eqtls_thyroid = rbind(Top_eqtls_thyroid,randomRows(genewise_min,1))
  }
  else
  {
    Top_eqtls_thyroid = rbind(Top_eqtls_thyroid,genewise_min)
  }
  
  rm(genewise)
  rm(genewise_min)
  message(k,"\r\n",appendLF=FALSE)
  flush.console()
  k=k+1
}