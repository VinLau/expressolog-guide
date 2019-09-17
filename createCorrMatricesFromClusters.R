#NB### SET YOUR WD HERE
setwd("~/mywd/")
#NB### 

#NB### INSTALL THESE PACKAGES IF NOT ALREADY INSTALLED
library("dplyr")
library("tidyverse")
library(data.table)
#NB###

#NB### NAME YOUR OUTPUT FILE NAME HERE (WARNING, THIS SCRIPT WILL DELETE A FILE NAMED THIS IN YOUR DIR)
flname <- "tom-ath-expressologCorMatrices.csv"
#NB### 
if (length(commandArgs(trailingOnly = TRUE)) >= 1){ flname <- commandArgs(trailingOnly = TRUE)[1][1] }

if (file.exists(flname)){ file.remove(flname) }

#NB### REPLACE THE NAME OF YOUR EXPRESSION FILE HERE
tomatoCSVDataTPM <- read.csv(file="20180920_Sl_ITAG3.2_Kallisto_quantile_norm_median_tpm_CompCT.csv",
                             header=TRUE, sep=",")
#NB### 

#NB### THESE LINES ARE FOR CLEANING TPM DATA, COMMENT THESE OUT IF YOU DONT NEED TO FILTER
zeroClndTomCSV <- tomatoCSVDataTPM[rowSums(tomatoCSVDataTPM[,-1]) > 0,]
fullyFltredTomCSV <- filter_all(zeroClndTomCSV, any_vars(. > 2))
#NB### 

#NB### THIS LINE REMOVES ROWNAMES, PARTICUARLY THE FIRST CELL
tomatoGenesExpr <- fullyFltredTomCSV %>% remove_rownames() %>% column_to_rownames(var="SlGeneID")
#NB### 

#NB### THIS LINE THEN REPLACES THE OLDER 'TISSUE' COLUMN NAMES WITH A UNIVERSAL ONE SO THAT IT IS EASIER FOR MERGING IN THE FUTURE
tomatoGenesExpr <- setnames(tomatoGenesExpr, old=c("Sl_COR","Sl_EN", "Sl_EP", "Sl_MCO", "Sl_MZ", "Sl_PH", "Sl_V", "Sl_QC", "Sl_35S"), new=c("COR", "EN", "EP", "MCO", "MZ", "PH", "V", "QC", "35S"))
#NB### 

#NB### IF YOU NEED TO DELETE A COLUMN (BECAUSE AN EQUILVALENT TISSUE DOES NOT EXIST IN THE OTHER EXPR DATA)
tomatoGenesExpr$QC = NULL
#NB### 

#NB### SAME LOGIC ABOVE APPLIED TO OTHER SPECIES' EXPR DATA
athCSVData <- read.csv(file="At_Mustroph_RMA_mean_unlogged_COMP_CT.csv", header=TRUE, sep=",")
athGenesExpr <- athCSVData %>% remove_rownames() %>% column_to_rownames(var="AtGeneID")
athGenesExpr <- setnames(athGenesExpr, old=c("At_RT_35S","At_RT_V", "At_RT_EN", "At_MZ", "At_EP", "At_PH", "At_MCO", "At_COR"), new=c("35S", "V", "EN", "MZ", "EP", "PH", "MCO", "COR"))
#NB### 

#NB### RENAME YOUR ORTHOMCL FILE HERE
omcl_read <- read.delim("jul_28_all_orthomcl.out", header = FALSE, sep = "\t")
omcl_read$V1 = NULL
#NB###

#NB### NAME YOUR CELL NAMES, WHAT YOU HAD NOTED ABOVE
cellNames <- c("COR", "EN", "MZ", "35S", "V", "EP", "PH", "MCO")
#NB### 

# Desired matrix for correlation matrix output
#       GeneA GeneA` GeneA```
# 35S  100  432    322
# EN   25    231    213
# I.e. the cell types are collated together to be the observations whilst Genes are variables

# Loop through every cluster
by(omcl_read, 1:nrow(omcl_read), function(row) {
  genes <- unlist(strsplit(as.character(row$V2), " "))[-1]
  
  # For now, delete rice-containing genes; if (grepl("IRGSP-1.0_protein_2019-06-26", gene)) { next }
  #NB### Below line deletes Rice containing genes in the oMCL cluster, add more of these if you need to delete more species
  genes <- genes[lapply(genes, function(x) length(grep("IRGSP-1.0_protein_2019-06-26", x, value=FALSE))) == 0]
  #NB### 
  
  # Create a empty matrix for correlation calculation
  geneCellMtx <- matrix(nrow=8,ncol=0)
  rownames(geneCellMtx) <- cellNames
  
  # loop through every gene in a cluster
  for(gene in genes){
    tempDF <- matrix()
    #NB###  Make sure to change the regex below to the species you want
    if (grepl("ITAG3.2_proteins", gene)) { 
      gene <- sub("\\(ITAG3.2_proteins\\)", "", gene)
      tempDF <- data.frame(t(tomatoGenesExpr[gene,]))
    }
    else if (grepl("TAIR10_pep_20101214_updated", gene)){
      gene <- sub("\\.\\d\\(TAIR10_pep_20101214_updated\\)", "", gene);
      tempDF <- data.frame(t(athGenesExpr[gene,]))
    }
    #NB### 
    
    if (anyNA(tempDF)) {next} # Check if we even have expr values for the gene
    
    # Check if column of expression values has already been added (e.g. A th expr for single gene for multiple transcripts)
    if (!colnames(tempDF)[1] %in% colnames(geneCellMtx)) { 
      geneCellMtx <- transform(merge(geneCellMtx, tempDF, by="row.names", all = TRUE ), row.names=Row.names, Row.names=NULL)
    }
    
    # print(tempDF)
  }
  # print(geneCellMtx)
  # print(cor(geneCellMtx))
  #NB### Set your correlation method here, default is PCC
  write.table(cor(geneCellMtx), sep = ",",file=flname, append=T, col.names=NA)
  #NB###
})

