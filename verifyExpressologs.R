#NB### Set your wd
setwd("~/mywd")
#NB###

#NB### Install packages if need be
library(data.table)
library("dplyr")
library("tidyverse")
library("ggplot2")
#NB###

#NB### You will need to change a lot of these variables, namely the file being read, the expression data filenames, and 'universal' tissue names
exprlogs <- read.csv(file="final.tsv", header=FALSE, sep="\t")

tomatoGenesExpr <- read.csv(file="20180920_Sl_ITAG3.2_Kallisto_quantile_norm_median_tpm_CompCT.csv",header=TRUE, sep=",")
tomatoGenesExpr <- tomatoGenesExpr %>% remove_rownames() %>% column_to_rownames(var="SlGeneID")
tomatoGenesExpr <- setnames(tomatoGenesExpr, old=c("Sl_COR","Sl_EN", "Sl_EP", "Sl_MCO", "Sl_MZ", "Sl_PH", "Sl_V", "Sl_QC", "Sl_35S"), new=c("COR", "EN", "EP", "MCO", "MZ", "PH", "V", "QC", "35S"))

athGenesExpr <- read.csv(file="At_Mustroph_RMA_mean_unlogged_COMP_CT.csv", header=TRUE, sep=",")
athGenesExpr <- athGenesExpr %>% remove_rownames() %>% column_to_rownames(var="AtGeneID")
athGenesExpr <- setnames(athGenesExpr, old=c("At_RT_35S","At_RT_V", "At_RT_EN", "At_MZ", "At_EP", "At_PH", "At_MCO", "At_COR"), new=c("35S", "V", "EN", "MZ", "EP", "PH", "MCO", "COR"))
#NB###

#NB### Use your fasta files available from your initial clustering to get their annotations, rename these
proteomeA <- readLines("TAIR10_pep_20101214_updated.fasta")
proteomeB <- readLines("ITAG3.2_proteins.fasta")
#NB###

#NB### This function uses a regexs (you should replace them with your own) to determine create a gene name from its fasta file. Note the extra processing for Arabidopsis is simply because the FA file for A th is split by '|'
makeGeneName <- function(geneName) {
  proteomeFasta <- proteomeA
  if (grepl("Solyc", geneName)){
    proteomeFasta <- proteomeB
    return (gsub("\\(.*\\)", "", proteomeFasta[grep(geneName, proteomeFasta)]))
  }
  return (paste(geneName, strsplit(proteomeFasta[grep(geneName,proteomeFasta)], "\\|")[[1]][3]))
}
#NB###

common_tissues <- intersect(colnames(athGenesExpr), colnames(tomatoGenesExpr))

trimmedExprlogs <- exprlogs[exprlogs$V1 %in% names(which(table(exprlogs$V1) >=6)), ] # trim exprlog table based on genes which have a minimum of X number of orthologs for manual inspection
geneDf <- trimmedExprlogs[trimmedExprlogs[, 1]==trimmedExprlogs[sample(nrow(trimmedExprlogs), 1), "V1"],] # grab a row filtered by a single gene in DF above randomly (sample(nrow))
geneDf <- geneDf[order(-geneDf$V3),] # reorder by expr cc
geneDf <- rbind(head(geneDf, 2), tail(geneDf, 2)) # get first 2 and alst 2 occurences

listOfGGplots <- list()
j <- 1
for( i in rownames(geneDf) ){ # loop through every exprlog matchup and plot them!
  print(geneDf)
  #NB### Use regex to determine which gene belogns to which
  speciesAGene <- toString(geneDf[i, "V1"])
  speciesBGene <- toString(geneDf[i, "V2"])
  if (grepl("Solyc", geneDf[i, "V1"])){
    speciesAGene <- toString(geneDf[i, "V2"])
    speciesBGene <- toString(geneDf[i, "V1"])
  }
  #NB###
  
  
  message('ath: ', speciesAGene, " tom: ", speciesBGene)
  exprDataMatchup <- as.data.frame(t(rbind(athGenesExpr[speciesAGene,common_tissues,], tomatoGenesExpr[speciesBGene ,common_tissues,])))
  print(exprDataMatchup)
  message(j)
  listOfGGplots[[j]] <- ggplot(exprDataMatchup, aes(x=!!ensym(speciesAGene), y=!!ensym(speciesBGene), label=rownames(exprDataMatchup))) +
    stat_smooth(method="lm", col = "red", se = FALSE, size=0.5) +
    geom_point() + geom_text(hjust=0.50, vjust=-0.55) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    labs(x = makeGeneName(speciesAGene), y = makeGeneName(speciesBGene))
  j <- j + 1
}

do.call(grid.arrange, listOfGGplots)