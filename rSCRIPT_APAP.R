# ----------------------------------------------------
## Script name: rScript
##
## Purpose of script: Analyze RNAseq DATA in 
## PAR1/PAR4 ECKO Mice with APAP Overdose
##
## Author: Rahul Rajala
## Date Created: 2023-09-25
##
## Copyright (c) Rahul Rajala, 2023
## Email: rahul-rajala@ouhsc.edu
##
## ---------------------------------------------------
## Notes:
##   
##
## ---------------------------------------------------

## set working directory for Mac and PC---------------
setwd("~/")
##-----------------------------------------------------

## Import Packages------------------------------------
library(dplyr)
library(patchwork)
library(VennDiagram)
library(BiocManager)
library(ggbiplot)
library(EnhancedVolcano)
library(plotly)

#####################Functions-----------------------
D_subset <- function(sam, BM)
{
  sam_ss <- subset(sam, baseMean > BM)
  sam_ss$PAdj <- p.adjust(sam_ss$PValue, method = "BH", n = nrow(sam_ss))
  return(sam_ss)
}

D_SIG <- function(sam1, Pcutoff)
{
  sam1_SIG <- subset(sam1, PAdj < Pcutoff)
}

heatprocess <- function(TSheet, GeneOI, Titles)
{
  raw_TSheet <- TSheet[, c(2, 14, 15, 16, 17, 18, 19)]
  SIG_TSheet <- GeneOI
  FR_K1 <- raw_TSheet[raw_TSheet$gene_name %in% SIG_TSheet, ]
  j <- FR_K1$gene_name
  FR_K1T <- as.data.frame(t(FR_K1[, -1]))
  colnames(FR_K1T) <- j
  FR_K1T$myfactor <- factor(row.names(FR_K1T))
  colnames(FR_K1T) <- j
  rownames(FR_K1T) <- Titles
  FR_K1T <- FR_K1T[, unlist(lapply(FR_K1T, is.numeric))]  # Remove non-numeric columns
  FR_K1T <- FR_K1T[, !apply(FR_K1T == 0, 2, all)]         # Remove non zero columns
  FR_K1T.pca <- prcomp(FR_K1T, center = TRUE, scale = TRUE)
  summary(FR_K1T.pca)
  str(FR_K1T.pca)
  #ggbiplot(FR_K1T.pca, labels = rownames(FR_K1T), var.axes = FALSE)
  nFR_K1T <- as.matrix(FR_K1T)
  heatmap(nFR_K1T, scale = "column", col = cm.colors(256))
}

PrePro_heatprocess <- function(Total, GeneOI, Titles)
{
  raw_TSheet <- Total
  SIG_TSheet <- GeneOI
  FR_K1 <- raw_TSheet[raw_TSheet$gene_name %in% SIG_TSheet, ]
  j <- FR_K1$gene_name
  FR_K1T <- as.data.frame(t(FR_K1[, -1]))
  colnames(FR_K1T) <- j
  FR_K1T$myfactor <- factor(row.names(FR_K1T))
  colnames(FR_K1T) <- j
  rownames(FR_K1T) <- Titles
  FR_K1T <- FR_K1T[, unlist(lapply(FR_K1T, is.numeric))]  # Remove non-numeric columns
  FR_K1T <- FR_K1T[, !apply(FR_K1T == 0, 2, all)]         # Remove non zero columns
  FR_K1T.pca <- prcomp(FR_K1T, center = TRUE, scale = TRUE)
  summary(FR_K1T.pca)
  str(FR_K1T.pca)
  #ggbiplot(FR_K1T.pca, labels = rownames(FR_K1T), var.axes = FALSE)
  nFR_K1T <- as.matrix(FR_K1T)
  heatmap(nFR_K1T, scale = "column", col = cm.colors(256))
}

#Feed this Function the datasheet, Genes of Interest, and Title of Samples and Variable AXES
PCAprocess <- function(TSheet, GeneOI, Titles, AXES)
{
  raw_TSheet <- TSheet[, c(2, 14, 15, 16, 17, 18, 19)]
  SIG_TSheet <- GeneOI
  FR_K1 <- raw_TSheet[raw_TSheet$gene_name %in% SIG_TSheet, ]
  j <- FR_K1$gene_name
  FR_K1T <- as.data.frame(t(FR_K1[, -1]))
  colnames(FR_K1T) <- j
  FR_K1T$myfactor <- factor(row.names(FR_K1T))
  colnames(FR_K1T) <- j
  rownames(FR_K1T) <- Titles
  FR_K1T <- FR_K1T[, unlist(lapply(FR_K1T, is.numeric))]  # Remove non-numeric columns
  FR_K1T <- FR_K1T[, !apply(FR_K1T == 0, 2, all)]         # Remove non zero columns
  FR_K1T.pca <- prcomp(FR_K1T, center = TRUE, scale = TRUE)
  summary(FR_K1T.pca)
  str(FR_K1T.pca)
  ggbiplot(FR_K1T.pca, labels = rownames(FR_K1T), var.axes = AXES)
}

PrePro_PCAprocess <- function(Total, GeneOI, Titles, AXES)
{
  FR_K1 <- Total
  j <- FR_K1[FR_K1$gene_name %in% GeneOI,]
  FR_K1T <- as.data.frame(t(FR_K1[,-1]))
  colnames(FR_K1T) <- j
  FR_K1T$myfactor <- factor(row.names(FR_K1T))
  colnames(FR_K1T) <- j
  rownames(FR_K1T) <- Titles
  FR_K1T <- FR_K1T[ , unlist(lapply(FR_K1T, is.numeric))]# Remove non-numeric columns
  FR_K1T <- FR_K1T[, !apply(FR_K1T == 0, 2, all)] # Remove non zero columns
  FR_K1T.pca <-prcomp(FR_K1T, center = TRUE, scale = TRUE)
  summary(FR_K1T.pca)
  str(FR_K1T.pca)
  ggbiplot(FR_K1T.pca, labels=rownames(FR_K1T), var.axes = FALSE)
}
##-----------------------------------------------------


#################Importing Data-DESEQ2-----------------
##Comparing Between Genotypes------------------------------------------
#APAP
dP1KABvP4KAB <- D_subset(read.csv("C:/Users/rahul/OneDrive/Desktop/HISAT_COMPARE/06-DE_lists/deseq2-Par1_ECko_APAP_Bead_vs_Par4_ECko_APAP_Bead_w_gene_names.tsv", sep = ""), 10)
dP1KABvP4KAB_SIG <- D_SIG(dP1KABvP4KAB, 0.05)

dP1KAIvP4KAI <- D_subset(read.csv("C:/Users/rahul/OneDrive/Desktop/HISAT_COMPARE/06-DE_lists/deseq2-Par1_ECko_APAP_Input_vs_Par4_ECko_APAP_Input_w_gene_names.tsv", sep = ""), 10)
dP1KAIvP4KAI_SIG <- D_SIG(dP1KAIvP4KAI, 0.05)

dP1HABvP4HAB <- D_subset(read.csv("C:/Users/rahul/OneDrive/Desktop/HISAT_COMPARE/06-DE_lists/deseq2-Par1_EChet_APAP_Bead_vs_Par4_EChet_APAP_Bead_w_gene_names.tsv", sep = ""), 10)
dP1HABvP4HAB_SIG <- D_SIG(dP1HABvP4HAB, 0.05)

dP1HAIvP4HAI <- D_subset(read.csv("C:/Users/rahul/OneDrive/Desktop/HISAT_COMPARE/06-DE_lists/deseq2-Par1_EChet_APAP_Input_vs_Par4_EChet_APAP_Input_w_gene_names.tsv", sep = ""), 10)
dP1HAIvP4HAI_SIG <- D_SIG(dP1HAIvP4HAI, 0.05)

#Saline
dP1KSBvP4KSB <- D_subset(read.csv("C:/Users/rahul/OneDrive/Desktop/HISAT_COMPARE/06-DE_lists/deseq2-Par1_ECko_Saline_Bead_vs_Par4_ECko_Saline_Bead_w_gene_names.tsv", sep = ""), 10)
dP1KSBvP4KSB_SIG <- D_SIG(dP1KSBvP4KSB, 0.05)

dP1KSIvP4KSI <- D_subset(read.csv("C:/Users/rahul/OneDrive/Desktop/HISAT_COMPARE/06-DE_lists/deseq2-Par1_ECko_Saline_Input_vs_Par4_ECko_Saline_Input_w_gene_names.tsv", sep = ""), 10)
dP1KSIvP4KSI_SIG <- D_SIG(dP1KSIvP4KSI, 0.05)

dP1HSBvP4HSB <- D_subset(read.csv("C:/Users/rahul/OneDrive/Desktop/HISAT_COMPARE/06-DE_lists/deseq2-Par1_EChet_Saline_Bead_vs_Par4_EChet_Saline_Bead_w_gene_names.tsv", sep = ""), 10)
dP1HSBvP4HSB_SIG <- D_SIG(dP1HSBvP4HSB, 0.05)

dP1HSIvP4HSI <- D_subset(read.csv("C:/Users/rahul/OneDrive/Desktop/HISAT_COMPARE/06-DE_lists/deseq2-Par1_EChet_Saline_Input_vs_Par4_EChet_Saline_Input_w_gene_names.tsv", sep = ""), 10)
dP1HSIvP4HSI_SIG <- D_SIG(dP1HSIvP4HSI, 0.05)

##Comparing Between KO------------------------------------------
#PAR1
dP1HSBvP1KSB <- D_subset(read.csv("C:/Users/rahul/OneDrive/Desktop/HISAT_COMPARE/06-DE_lists/deseq2-Par1_EChet_Saline_Bead_vs_Par1_ECko_Saline_Bead_w_gene_names.tsv", sep = ""), 10)
dP1HSBvP1KSB_SIG <- D_SIG(dP1HSBvP1KSB, 0.05)

dP1HABvP1KAB <- D_subset(read.csv("C:/Users/rahul/OneDrive/Desktop/HISAT_COMPARE/06-DE_lists/deseq2-Par1_EChet_APAP_Bead_vs_Par1_ECko_APAP_Bead_w_gene_names.tsv", sep = ""), 10)
dP1HABvP1KAB_SIG <- D_SIG(dP1HABvP1KAB, 0.05)
heatprocess(dP1HABvP1KAB, dP1HABvP1KAB_SIG[, "gene_name"], c("Par1 EChet APAP Bead 1","Par1 EChet APAP Bead 2","Par1 EChet APAP Bead 3","Par1 ECko APAP Bead 1","Par1 ECko APAP Bead 2","Par1 ECko APAP Bead 3"))

dP1HSIvP1KSI <- D_subset(read.csv("C:/Users/rahul/OneDrive/Desktop/HISAT_COMPARE/06-DE_lists/deseq2-Par1_EChet_Saline_Input_vs_Par1_ECko_Saline_Input_w_gene_names.tsv", sep = ""), 10)
dP1HSIvP1KSI_SIG <- D_SIG(dP1HSIvP1KSI, 0.05)

dP1HAIvP1KAI <- D_subset(read.csv("C:/Users/rahul/OneDrive/Desktop/HISAT_COMPARE/06-DE_lists/deseq2-Par1_EChet_APAP_Input_vs_Par1_ECko_APAP_Input_w_gene_names.tsv", sep = ""), 10)
dP1HAIvP1KAI_SIG <- D_SIG(dP1HAIvP1KAI, 0.05)

#PAR4
dP4HSBvP4KSB <- D_subset(read.csv("C:/Users/rahul/OneDrive/Desktop/HISAT_COMPARE/06-DE_lists/deseq2-Par4_EChet_Saline_Bead_vs_Par4_ECko_Saline_Bead_w_gene_names.tsv", sep = ""), 10)
dP4HSBvP4KSB_SIG <- D_SIG(dP4HSBvP4KSB, 0.05)

dP4HABvP4KAB <- D_subset(read.csv("C:/Users/rahul/OneDrive/Desktop/HISAT_COMPARE/06-DE_lists/deseq2-Par4_EChet_APAP_Bead_vs_Par4_ECko_APAP_Bead_w_gene_names.tsv", sep = ""), 10)
dP4HABvP4KAB_SIG <- D_SIG(dP4HABvP4KAB, 0.05)

dP4HSIvP4KSI <- D_subset(read.csv("C:/Users/rahul/OneDrive/Desktop/HISAT_COMPARE/06-DE_lists/deseq2-Par4_EChet_Saline_Input_vs_Par4_ECko_Saline_Input_w_gene_names.tsv", sep = ""), 10)
dP4HSIvP4KSI_SIG <- D_SIG(dP4HSIvP4KSI, 0.05)

dP4HAIvP4KAI <- D_subset(read.csv("C:/Users/rahul/OneDrive/Desktop/HISAT_COMPARE/06-DE_lists/deseq2-Par4_EChet_APAP_Input_vs_Par4_ECko_APAP_Input_w_gene_names.tsv", sep = ""), 10)
dP4HAIvP4KAI_SIG <- D_SIG(dP4HAIvP4KAI, 0.05)

##Comparing Between Treatment------------------------------------------
#PAR1
dP1HSBvP1HAB <- D_subset(read.csv("C:/Users/rahul/OneDrive/Desktop/HISAT_COMPARE/06-DE_lists/deseq2-Par1_EChet_Saline_Bead_vs_Par1_EChet_APAP_Bead_w_gene_names.tsv", sep = ""), 10)
dP1HSBvP1HAB_SIG <- D_SIG(dP1HSBvP1HAB, 0.05)

dP1HSIvP1HAI <- D_subset(read.csv("C:/Users/rahul/OneDrive/Desktop/HISAT_COMPARE/06-DE_lists/deseq2-Par1_EChet_Saline_Input_vs_Par1_EChet_APAP_Input_w_gene_names.tsv", sep = ""), 10)
dP1HSIvP1HAI_SIG <- D_SIG(dP1HSIvP1HAI, 0.05)

dP1KSBvP1KAB <- D_subset(read.csv("C:/Users/rahul/OneDrive/Desktop/HISAT_COMPARE/06-DE_lists/deseq2-Par1_ECko_Saline_Bead_vs_Par1_ECko_APAP_Bead_w_gene_names.tsv", sep = ""), 10)
dP1KSBvP1KAB_SIG <- D_SIG(dP1KSBvP1KAB, 0.05)

dP1KSIvP1KAI <- D_subset(read.csv("C:/Users/rahul/OneDrive/Desktop/HISAT_COMPARE/06-DE_lists/deseq2-Par1_ECko_Saline_Input_vs_Par1_ECko_APAP_Input_w_gene_names.tsv", sep = ""), 10)
dP1KSIvP1KAI_SIG <- D_SIG(dP1KSIvP1KAI, 0.05)

#PAR4
dP4HSBvP4HAB <- D_subset(read.csv("C:/Users/rahul/OneDrive/Desktop/HISAT_COMPARE/06-DE_lists/deseq2-Par4_EChet_Saline_Bead_vs_Par4_EChet_APAP_Bead_w_gene_names.tsv", sep = ""), 10)
dP4HSBvP4HAB_SIG <- D_SIG(dP4HSBvP4HAB, 0.05)

dP4HSIvP4HAI <- D_subset(read.csv("C:/Users/rahul/OneDrive/Desktop/HISAT_COMPARE/06-DE_lists/deseq2-Par4_EChet_Saline_Input_vs_Par4_EChet_APAP_Input_w_gene_names.tsv", sep = ""), 10)
dP4HSIvP4HAI_SIG <- D_SIG(dP4HSIvP4HAI, 0.05)

dP4KSBvP4KAB <- D_subset(read.csv("C:/Users/rahul/OneDrive/Desktop/HISAT_COMPARE/06-DE_lists/deseq2-Par4_ECko_Saline_Bead_vs_Par4_ECko_APAP_Bead_w_gene_names.tsv", sep = ""), 10)
dP4KSBvP4KAB_SIG <- D_SIG(dP4KSBvP4KAB, 0.05)

dP4KSIvP4KAI <- D_subset(read.csv("C:/Users/rahul/OneDrive/Desktop/HISAT_COMPARE/06-DE_lists/deseq2-Par4_ECko_Saline_Input_vs_Par4_ECko_APAP_Input_w_gene_names.tsv", sep = ""), 10)
dP4KSIvP4KAI_SIG <- D_SIG(dP4KSIvP4KAI, 0.05)

##Comparing Between Fraction------------------------------------------
#PAR1
dP1HSIvP1HSB <- D_subset(read.csv("C:/Users/rahul/OneDrive/Desktop/HISAT_COMPARE/06-DE_lists/deseq2-Par1_EChet_Saline_Input_vs_Par1_EChet_Saline_Bead_w_gene_names.tsv", sep = ""), 10)
dP1HSIvP1HSB_SIG <- D_SIG(dP1HSIvP1HSB, 0.05)

dP1KSIvP1KSB <- D_subset(read.csv("C:/Users/rahul/OneDrive/Desktop/HISAT_COMPARE/06-DE_lists/deseq2-Par1_ECko_Saline_Input_vs_Par1_ECko_Saline_Bead_w_gene_names.tsv", sep = ""), 10)
dP1KSIvP1KSB_SIG <- D_SIG(dP1KSIvP1KSB, 0.05)

dP1HAIvP1HAB <- D_subset(read.csv("C:/Users/rahul/OneDrive/Desktop/HISAT_COMPARE/06-DE_lists/deseq2-Par1_EChet_APAP_Input_vs_Par1_EChet_APAP_Bead_w_gene_names.tsv", sep = ""), 10)
dP1HAIvP1HAB_SIG <- D_SIG(dP1HAIvP1HAB, 0.05)

dP1KAIvP1KAB <- D_subset(read.csv("C:/Users/rahul/OneDrive/Desktop/HISAT_COMPARE/06-DE_lists/deseq2-Par1_ECko_APAP_Input_vs_Par1_ECko_APAP_Bead_w_gene_names.tsv", sep = ""), 10)
dP1KAIvP1KAB_SIG <- D_SIG(dP1KAIvP1KAB, 0.05)

#PAR4
dP4HSIvP4HSB <- D_subset(read.csv("C:/Users/rahul/OneDrive/Desktop/HISAT_COMPARE/06-DE_lists/deseq2-Par4_EChet_Saline_Input_vs_Par4_EChet_Saline_Bead_w_gene_names.tsv", sep = ""), 10)
dP4HSIvP4HSB_SIG <- D_SIG(dP4HSIvP4HSB, 0.05)

dP4KSIvP4KSB <- D_subset(read.csv("C:/Users/rahul/OneDrive/Desktop/HISAT_COMPARE/06-DE_lists/deseq2-Par4_ECko_Saline_Input_vs_Par4_ECko_Saline_Bead_w_gene_names.tsv", sep = ""), 10)
dP4KSIvP4KSB_SIG <- D_SIG(dP4KSIvP4KSB, 0.05)

dP4HAIvP4HAB <- D_subset(read.csv("C:/Users/rahul/OneDrive/Desktop/HISAT_COMPARE/06-DE_lists/deseq2-Par4_EChet_APAP_Input_vs_Par4_EChet_APAP_Bead_w_gene_names.tsv", sep = ""), 10)
dP4HAIvP4HAB_SIG <- D_SIG(dP4HAIvP4HAB, 0.05)

dP4KAIvP4KAB <- D_subset(read.csv("C:/Users/rahul/OneDrive/Desktop/HISAT_COMPARE/06-DE_lists/deseq2-Par4_ECko_APAP_Input_vs_Par4_ECko_APAP_Bead_w_gene_names.tsv", sep = ""), 10)
dP4KAIvP4KAB_SIG <- D_SIG(dP4KAIvP4KAB, 0.05)

dP1HSIvP1KAI <- D_subset(read.csv("C:/Users/rahul/OneDrive/Desktop/HISAT_COMPARE/06-DE_lists/deseq2-Par1_EChet_Saline_Input_vs_Par1_ECko_APAP_Input_w_gene_names.tsv", sep = ""), 10)
dP1HSIvP1KAI_SIG <- D_SIG(dP1HSIvP1KAI, 0.05)
write.csv(x = dP1HSIvP1KAI, file = "Par1_EChet_Saline_Input_vs_Par1_ECko_APAP_Input.csv")

dP1HSIvP1KAI <- D_subset(read.csv("C:/Users/rahul/OneDrive/Desktop/HISAT_COMPARE/06-DE_lists/deseq2-Par1_EChet_Saline_Input_vs_Par1_ECko_APAP_Input_w_gene_names.tsv", sep = ""), 10)
dP1HSIvP1KAI_SIG <- D_SIG(dP1HSIvP1KAI, 0.05)
write.csv(x = dP1HSIvP1KAI, file = "Par1_EChet_Saline_Input_vs_Par1_ECko_APAP_Input.csv")

##-----------------------------------------------------


###############################Volcano Plots-----------
##Comparing Between Genotypes------------------------------------------
#APAP
EnhancedVolcano(dP1KABvP4KAB, lab = dP1KABvP4KAB$gene_name , x = 'log2FoldChange', y = 'PAdj', FCcutoff = 0.5, pCutoff = 0.05, title = "Par1 ECko APAP Bead vs Par4 ECko APAP Bead")
EnhancedVolcano(dP1KAIvP4KAI, lab = dP1KAIvP4KAI$gene_name , x = 'log2FoldChange', y = 'PAdj', FCcutoff = 0.5, pCutoff = 0.05, title = "Par1 ECko APAP Input vs Par4 ECko APAP Input")
EnhancedVolcano(dP1HABvP4HAB, lab = dP1HABvP4HAB$gene_name , x = 'log2FoldChange', y = 'PAdj', FCcutoff = 0.5, pCutoff = 0.05, title = "Par1 EChet APAP Bead vs Par4 EChet APAP Bead")
EnhancedVolcano(dP1HAIvP4HAI, lab = dP1HAIvP4HAI$gene_name , x = 'log2FoldChange', y = 'PAdj', FCcutoff = 0.5, pCutoff = 0.05, title = "Par1 EChet APAP Input vs Par4 EChet APAP Input")

#Saline
EnhancedVolcano(dP1HSBvP1KSB, lab = dP1HSBvP1KSB$gene_name , x = 'log2FoldChange', y = 'PAdj', FCcutoff = 0.5, pCutoff = 0.05, title = "Par1 ECko Saline Bead vs Par4 ECko Saline Bead")
EnhancedVolcano(dP1KSIvP4KSI, lab = dP1KSIvP4KSI$gene_name , x = 'log2FoldChange', y = 'PAdj', FCcutoff = 0.5, pCutoff = 0.05, title = "Par1 ECko Saline Input vs Par4 ECko Saline Input")
EnhancedVolcano(dP1HSBvP4HSB, lab = dP1HSBvP4HSB$gene_name , x = 'log2FoldChange', y = 'PAdj', FCcutoff = 0.5, pCutoff = 0.05, title = "Par1 EChet Saline Bead vs Par4 EChet Saline Bead")
EnhancedVolcano(dP1HSBvP4HSB, lab = dP1HSBvP4HSB$gene_name , x = 'log2FoldChange', y = 'PAdj', FCcutoff = 0.5, pCutoff = 0.05, title = "Par1 EChet Saline Input vs Par4 EChet Saline Input")

##Comparing Between KO------------------------------------------
#PAR1
EnhancedVolcano(dP1HSBvP1KSB, lab = dP1HSBvP1KSB$gene_name , x = 'log2FoldChange', y = 'PAdj', FCcutoff = 0.5, pCutoff = 0.05, title = "Par1 EChet Saline Bead vs Par1 ECko Saline Bead")
EnhancedVolcano(dP1HABvP1KAB, lab = dP1HABvP1KAB$gene_name , x = 'log2FoldChange', y = 'PAdj', FCcutoff = 0.5, pCutoff = 0.05, title = "Par1 EChet APAP Bead vs Par1 ECko APAP Bead")
EnhancedVolcano(dP1HSIvP1KSI, lab = dP1HSIvP1KSI$gene_name , x = 'log2FoldChange', y = 'PAdj', FCcutoff = 0.5, pCutoff = 0.05, title = "Par1 EChet Saline Input vs Par1 ECko Saline Input")
EnhancedVolcano(dP1HAIvP1KAI, lab = dP1HAIvP1KAI$gene_name , x = 'log2FoldChange', y = 'PAdj', FCcutoff = 0.5, pCutoff = 0.05, title = "Par1 EChet APAP Input vs Par1 ECko APAP Input")

#PAR4
EnhancedVolcano(dP4HSBvP4KSB, lab = dP4HSBvP4KSB$gene_name , x = 'log2FoldChange', y = 'PAdj', FCcutoff = 0.5, pCutoff = 0.05, title = "Par4 EChet Saline Bead vs Par4 ECko Saline Bead")
EnhancedVolcano(dP4HABvP4KAB, lab = dP4HABvP4KAB$gene_name , x = 'log2FoldChange', y = 'PAdj', FCcutoff = 0.5, pCutoff = 0.05, title = "Par4 EChet APAP Bead vs Par4 ECko APAP Bead")
EnhancedVolcano(dP4HSIvP4KSI, lab = dP4HSIvP4KSI$gene_name , x = 'log2FoldChange', y = 'PAdj', FCcutoff = 0.5, pCutoff = 0.05, title = "Par4 EChet Saline Input vs Par4 ECko Saline Input")
EnhancedVolcano(dP4HAIvP4KAI, lab = dP4HAIvP4KAI$gene_name , x = 'log2FoldChange', y = 'PAdj', FCcutoff = 0.5, pCutoff = 0.05, title = "Par4 EChet APAP Input vs Par4 ECko APAP Input")

##Comparing Between Treatment------------------------------------------
#PAR1
EnhancedVolcano(dP1HSBvP1HAB, lab = dP1HSBvP1HAB$gene_name , x = 'log2FoldChange', y = 'PAdj', FCcutoff = 0.5, pCutoff = 0.05, title = "Par1 EChet Saline Bead vs Par1 EChet APAP Bead")
EnhancedVolcano(dP1HSIvP1HAI, lab = dP1HSIvP1HAI$gene_name , x = 'log2FoldChange', y = 'PAdj', FCcutoff = 0.5, pCutoff = 0.05, title = "Par1 EChet Saline Input vs Par1 EChet APAP Input")
EnhancedVolcano(dP1KSBvP1KAB, lab = dP1KSBvP1KAB$gene_name , x = 'log2FoldChange', y = 'PAdj', FCcutoff = 0.5, pCutoff = 0.05, title = "Par1 ECko Saline Bead vs Par1 ECko APAP Bead")
EnhancedVolcano(dP1KSIvP1KAI, lab = dP1KSIvP1KAI$gene_name , x = 'log2FoldChange', y = 'PAdj', FCcutoff = 0.5, pCutoff = 0.05, title = "Par1 ECko Saline Input vs Par1 ECko APAP Input")

#PAR4
EnhancedVolcano(dP4HSBvP4HAB, lab = dP4HSBvP4HAB$gene_name , x = 'log2FoldChange', y = 'PAdj', FCcutoff = 0.5, pCutoff = 0.05, title = "Par4 EChet Saline Bead vs Par4 EChet APAP Bead")
EnhancedVolcano(dP4HSIvP4HAI, lab = dP4HSIvP4HAI$gene_name , x = 'log2FoldChange', y = 'PAdj', FCcutoff = 0.5, pCutoff = 0.05, title = "Par4 EChet Saline Input vs Par4 EChet APAP Input")
EnhancedVolcano(dP4KSBvP4KAB, lab = dP4KSBvP4KAB$gene_name , x = 'log2FoldChange', y = 'PAdj', FCcutoff = 0.5, pCutoff = 0.05, title = "Par4 ECko Saline Bead vs Par4 ECko APAP Bead")
EnhancedVolcano(dP4KSIvP4KAI, lab = dP4KSIvP4KAI$gene_name , x = 'log2FoldChange', y = 'PAdj', FCcutoff = 0.5, pCutoff = 0.05, title = "Par4 ECko Saline Input vs Par4 ECko APAP Input")

##Comparing Between Fraction------------------------------------------
#PAR1
EnhancedVolcano(dP1HSIvP1HSB, lab = dP1HSIvP1HSB$gene_name , x = 'log2FoldChange', y = 'PAdj', FCcutoff = 0.5, pCutoff = 0.05, title = "Par1 EChet Saline Input vs Par1 EChet Saline Bead")
EnhancedVolcano(dP1KSIvP1KSB, lab = dP1KSIvP1KSB$gene_name , x = 'log2FoldChange', y = 'PAdj', FCcutoff = 0.5, pCutoff = 0.05, title = "Par1 ECko Saline Input vs Par1 ECko Saline Bead")
EnhancedVolcano(dP1HAIvP1HAB, lab = dP1HAIvP1HAB$gene_name , x = 'log2FoldChange', y = 'PAdj', FCcutoff = 0.5, pCutoff = 0.05, title = "Par1 EChet APAP Input vs Par1 EChet APAP Bead")
EnhancedVolcano(dP1KAIvP1KAB, lab = dP1KAIvP1KAB$gene_name , x = 'log2FoldChange', y = 'PAdj', FCcutoff = 0.5, pCutoff = 0.05, title = "Par1 ECko APAP Input vs Par1 ECko APAP Bead")

#PAR4
EnhancedVolcano(dP4HSIvP4HSB, lab = dP4HSIvP4HSB$gene_name , x = 'log2FoldChange', y = 'PAdj', FCcutoff = 0.5, pCutoff = 0.05, title = "Par4 EChet Saline Input vs Par4 EChet Saline Bead")
EnhancedVolcano(dP4KSIvP4KSB, lab = dP4KSIvP4KSB$gene_name , x = 'log2FoldChange', y = 'PAdj', FCcutoff = 0.5, pCutoff = 0.05, title = "Par4 ECko Saline Input vs Par4 ECko Saline Bead")
EnhancedVolcano(dP4HAIvP4HAB, lab = dP4HAIvP4HAB$gene_name , x = 'log2FoldChange', y = 'PAdj', FCcutoff = 0.5, pCutoff = 0.05, title = "Par4 EChet APAP Input vs Par4 EChet APAP Bead")
EnhancedVolcano(dP4KAIvP4KAB, lab = dP4KAIvP4KAB$gene_name , x = 'log2FoldChange', y = 'PAdj', FCcutoff = 0.5, pCutoff = 0.05, title = "Par4 ECko APAP Input vs Par4 ECko APAP Bead")
##------------------------------------------------------

##------------------------------------------------------
##Heatmaps----------------------------------------------
##Comparing Between Genotypes------------------------------------------
#APAP
heatprocess(dP1KABvP4KAB, dP1KABvP4KAB_SIG[, "gene_name"], c("Par1 ECko APAP Bead 1","Par1 ECko APAP Bead 2","Par1 ECko APAP Bead 3","Par4 ECko APAP Bead 1","Par4 ECko APAP Bead 2","Par4 ECko APAP Bead 3"))
heatprocess(dP1KAIvP4KAI, dP1KAIvP4KAI_SIG[, "gene_name"], c("Par1 ECko APAP Input 1","Par1 ECko APAP Input 2","Par1 ECko APAP Input 3","Par4 ECko APAP Input 1","Par4 ECko APAP Input 2","Par4 ECko APAP Input 3"))
heatprocess(dP1HABvP4HAB, dP1HABvP4HAB_SIG[, "gene_name"], c("Par1 EChet APAP Bead 1","Par1 EChet APAP Bead 2","Par1 EChet APAP Bead 3","Par4 EChet APAP Bead 1","Par4 EChet APAP Bead 2","Par4 EChet APAP Bead 3"))
heatprocess(dP1HAIvP4HAI, dP1HAIvP4HAI_SIG[, "gene_name"], c("Par1 EChet APAP Input 1","Par1 EChet APAP Input 2","Par1 EChet APAP Input 3","Par4 EChet APAP Input 1","Par4 EChet APAP Input 2","Par4 EChet APAP Input 3"))

#Saline
heatprocess(dP1KSBvP4KSB, dP1KSBvP4KSB_SIG[, "gene_name"], c("Par1 ECko Saline Bead 1","Par1 ECko Saline Bead 2","Par1 ECko Saline Bead 3","Par4 ECko Saline 1","Par4 ECko Saline 2","Par4 ECko Saline 3"))
heatprocess(dP1KSIvP4KSI, dP1KSIvP4KSI_SIG[, "gene_name"], c("Par1 ECko Saline Input 1","Par1 ECko Saline Input 2","Par1 ECko Saline Input 3","Par4 ECko Saline Input 1","Par4 ECko Saline Input 2","Par4 ECko Saline Input 3"))
heatprocess(dP1HSBvP4HSB, dP1HSBvP4HSB_SIG[, "gene_name"], c("Par1 EChet Saline Bead 1","Par1 EChet Saline Bead 2","Par1 EChet Saline Bead 3","Par4 EChet Saline Bead 1","Par4 EChet Saline Bead 2","Par4 EChet Saline Bead 3"))
heatprocess(dP1HSIvP4HSI, dP1HSIvP4HSI_SIG[, "gene_name"], c("Par1 EChet Saline Input 1","Par1 EChet Saline Input 2","Par1 EChet Saline Input 3","Par4 EChet Saline Input 1","Par4 EChet Saline Input 2","Par4 EChet Saline Input 3"))

##Comparing Between KO------------------------------------------
#PAR1
heatprocess(dP1HSBvP1KSB, dP1HSBvP1KSB_SIG[, "gene_name"], c("Par1 EChet Saline Bead 1","Par1 EChet Saline Bead 2","Par1 EChet Saline Bead 3","Par1 ECko Saline Bead 1","Par1 ECko Saline Bead 2","Par1 ECko Saline Bead 3"))
heatprocess(dP1HABvP1KAB, dP1HABvP1KAB_SIG[, "gene_name"], c("Par1 EChet APAP Bead 1","Par1 EChet APAP Bead 2","Par1 EChet APAP Bead 3","Par1 ECko APAP Bead 1","Par1 ECko APAP Bead 2","Par1 ECko APAP Bead 3"))
heatprocess(dP1HSIvP1KSI, dP1HSIvP1KSI_SIG[, "gene_name"], c("Par1 EChet Saline Input 1","Par1 EChet Saline Input 2","Par1 EChet Saline Input 3","Par1 ECko Saline Input 1","Par1 ECko Saline Input 2","Par1 ECko Saline Input 3"))
heatprocess(dP1HAIvP1KAI, dP1HAIvP1KAI_SIG[, "gene_name"], c("Par1 EChet APAP Input 1","Par1 EChet APAP Input 2","Par1 EChet APAP Input 3","Par1 ECko APAP Input 1","Par1 ECko APAP Input 2","Par1 ECko APAP Input 3"))

#PAR4
heatprocess(dP4HSBvP4KSB, dP4HSBvP4KSB_SIG[, "gene_name"], c("Par4 EChet Saline Bead 1","Par4 EChet Saline Bead 2","Par4 EChet Saline Bead 3","Par4 ECko Saline Bead 1","Par4 ECko Saline Bead 2","Par4 ECko Saline Bead 3"))
heatprocess(dP4HABvP4KAB, dP4HABvP4KAB_SIG[, "gene_name"], c("Par4 EChet APAP Bead 1","Par4 EChet APAP Bead 2","Par4 EChet APAP Bead 3","Par4 ECko APAP Bead 1","Par4 ECko APAP Bead 2","Par4 ECko APAP Bead 3"))
heatprocess(dP4HSIvP4KSI, dP4HSIvP4KSI_SIG[, "gene_name"], c("Par4 EChet Saline Input 1","Par4 EChet Saline Input 2","Par4 EChet Saline Input 3","Par4 ECko Saline Input 1","Par4 ECko Saline Input 2","Par4 ECko Saline Input 3"))
heatprocess(dP4HAIvP4KAI, dP4HAIvP4KAI_SIG[, "gene_name"], c("Par4 EChet APAP Input 1","Par4 EChet APAP Input 2","Par4 EChet APAP Input 3","Par4 ECko APAP Input 1","Par4 ECko APAP Input 2","Par4 ECko APAP Input 3"))

##Comparing Between Treatment------------------------------------------
#PAR1
heatprocess(dP1HSBvP1HAB, dP1HSBvP1HAB_SIG[, "gene_name"], c("Par1 EChet Saline Bead 1","Par1 EChet Saline Bead 2","Par1 EChet Saline Bead 3","Par1 EChet APAP Bead 1","Par1 EChet APAP Bead 2","Par1 EChet APAP Bead 3"))
heatprocess(dP1HSIvP1HAI, dP1HSIvP1HAI_SIG[, "gene_name"], c("Par1 EChet Saline Input 1","Par1 EChet Saline Input 2","Par1 EChet Saline Input 3","Par1 EChet APAP Input 1","Par1 EChet APAP Input 2","Par1 EChet APAP Input 3"))
heatprocess(dP1KSBvP1KAB, dP1KSBvP1KAB_SIG[, "gene_name"], c("Par1 ECko Saline Bead 1","Par1 ECko Saline Bead 2","Par1 ECko Saline Bead 3","Par1 ECko APAP Bead 1","Par1 ECko APAP Bead 2","Par1 ECko APAP Bead 3"))
heatprocess(dP1KSIvP1KAI, dP1KSIvP1KAI_SIG[, "gene_name"], c("Par1 ECko Saline Input 1","Par1 ECko Saline Input 2","Par1 ECko Saline Input 3","Par1 ECko APAP Input 1","Par1 ECko APAP Input 2","Par1 ECko APAP Input 3"))

#PAR4
heatprocess(dP4HSBvP4HAB, dP4HSBvP4HAB_SIG[, "gene_name"], c("Par4 EChet Saline Bead 1","Par4 EChet Saline Bead 2","Par4 EChet Saline Bead 3","Par4 EChet APAP Bead 1","Par4 EChet APAP Bead 2","Par4 EChet APAP Bead 3"))
heatprocess(dP4HSIvP4HAI, dP4HSIvP4HAI_SIG[, "gene_name"], c("Par4 EChet Saline Input 1","Par4 EChet Saline Input 2","Par4 EChet Saline Input 3","Par4 EChet APAP Input 1","Par4 EChet APAP Input 2","Par4 EChet APAP Input 3"))
heatprocess(dP4KSBvP4KAB, dP4KSBvP4KAB_SIG[, "gene_name"], c("Par4 ECko Saline Bead 1","Par4 ECko Saline Bead 2","Par4 ECko Saline Bead 3","Par4 ECko APAP Bead 1","Par4 ECko APAP Bead 2","Par4 ECko APAP Bead 3"))
heatprocess(dP4KSIvP4KAI, dP4KSIvP4KAI_SIG[, "gene_name"], c("Par4 ECko Saline Input 1","Par4 ECko Saline Input 2","Par4 ECko Saline Input 3","Par4 ECko APAP Input 1","Par4 ECko APAP Input 2","Par4 ECko APAP Input 3"))

##Comparing Between Fraction------------------------------------------
#PAR1
heatprocess(dP1HSIvP1HSB, dP1HSIvP1HSB_SIG[, "gene_name"], c("Par1 EChet Saline Input 1","Par1 EChet Saline Input 2","Par1 EChet Saline Input 3","Par1 EChet Saline Bead 1","Par1 EChet Saline Bead 2","Par1 EChet Saline Bead 3"))
heatprocess(dP1KSIvP1KSB, dP1KSIvP1KSB_SIG[, "gene_name"], c("Par1 ECko Saline Input 1","Par1 ECko Saline Input 2","Par1 ECko Saline Input 3","Par1 ECko Saline Bead 1","Par1 ECko Saline Bead 2","Par1 ECko Saline Bead 3"))
heatprocess(dP1HAIvP1HAB, dP1HAIvP1HAB_SIG[, "gene_name"], c("Par1 EChet APAP Input 1","Par1 EChet APAP Input 2","Par1 EChet APAP Input 3","Par1 EChet APAP Bead 1","Par1 EChet APAP Bead 2","Par1 EChet APAP Bead 3"))
heatprocess(dP1KAIvP1KAB, dP1KAIvP1KAB_SIG[, "gene_name"], c("Par1 ECko APAP Input 1","Par1 ECko APAP Input 2","Par1 ECko APAP Input 3","Par1 ECko APAP Bead 1","Par1 ECko APAP Bead 2","Par1 ECko APAP Bead 3"))

#PAR4
heatprocess(dP4HSIvP4HSB, dP4HSIvP4HSB_SIG[, "gene_name"], c("Par4 EChet Saline Input 1","Par4 EChet Saline Input 2","Par4 EChet Saline Input 3","Par4 EChet Saline Bead 1","Par4 EChet Saline Bead 2","Par4 EChet Saline Bead 3"))
heatprocess(dP4KSIvP4KSB, dP4KSIvP4KSB_SIG[, "gene_name"], c("Par4 ECko Saline Input 1","Par4 ECko Saline Input 2","Par4 ECko Saline Input 3","Par4 ECko Saline Bead 1","Par4 ECko Saline Bead 2","Par4 ECko Saline Bead 3"))
heatprocess(dP4HAIvP4HAB, dP4HAIvP4HAB_SIG[, "gene_name"], c("Par4 EChet APAP Input 1","Par4 EChet APAP Input 2","Par4 EChet APAP Input 3","Par4 EChet APAP Bead 1","Par4 EChet APAP Bead 2","Par4 EChet APAP Bead 3"))
heatprocess(dP4KAIvP4KAB, dP4KAIvP4KAB_SIG[, "gene_name"], c("Par4 ECko APAP Input 1","Par4 ECko APAP Input 2","Par4 ECko APAP Input 3","Par4 ECko APAP Bead 1","Par4 ECko APAP Bead 2","Par4 ECko Bead Input 3"))

##-----------------------------------------------------


##-----------------------------------------------------
##Total Datasheet #PAR1-HET-------------------
Total_Input_P1 <- read.csv("C:/Users/rahul/OneDrive/Desktop/HISAT_COMPARE/06-DE_lists/deseq2-Par1_EChet_Saline_Input_vs_Par1_EChet_APAP_Input_w_gene_names.tsv", sep = "")
Total_Bead_P1 <- read.csv("C:/Users/rahul/OneDrive/Desktop/HISAT_COMPARE/06-DE_lists/deseq2-Par1_EChet_Saline_Bead_vs_Par1_EChet_APAP_Bead_w_gene_names.tsv", sep = "")

TotalPAR1_HET <- cbind(Total_Input_P1[, c(2, 14, 15, 16, 17, 18, 19)], Total_Bead_P1[, c(2, 14, 15, 16, 17, 18, 19)])
TotalPAR1_HET <- TotalPAR1_HET[ , -c(8)]

##Total Datasheet #PAR4-HET------------------------------
Total_Input_P4 <- read.csv("C:/Users/rahul/OneDrive/Desktop/HISAT_COMPARE/06-DE_lists/deseq2-Par4_EChet_Saline_Input_vs_Par4_EChet_APAP_Input_w_gene_names.tsv", sep = "")
Total_Bead_P4 <- read.csv("C:/Users/rahul/OneDrive/Desktop/HISAT_COMPARE/06-DE_lists/deseq2-Par4_EChet_Saline_Bead_vs_Par4_EChet_APAP_Bead_w_gene_names.tsv", sep = "")

TotalPAR4_HET <- cbind(Total_Input_P4[, c(2, 14, 15, 16, 17, 18, 19)], Total_Bead_P4[, c(2, 14, 15, 16, 17, 18, 19)])
TotalPAR4_HET <- TotalPAR4_HET[ , -c(8)]

Total <- cbind(TotalPAR1_HET, TotalPAR4_HET)
Total <- Total[ , -c(14)]

##Total Datasheet #PAR1-KO------------------
Total_Input_P1_K <- read.csv("C:/Users/rahul/OneDrive/Desktop/HISAT_COMPARE/06-DE_lists/deseq2-Par1_ECko_Saline_Input_vs_Par1_ECko_APAP_Input_w_gene_names.tsv", sep = "")
Total_Bead_P1_K <- read.csv("C:/Users/rahul/OneDrive/Desktop/HISAT_COMPARE/06-DE_lists/deseq2-Par1_ECko_Saline_Bead_vs_Par1_ECko_APAP_Bead_w_gene_names.tsv", sep = "")

TotalPAR1_KO <- cbind(Total_Input_P1_K[, c(2, 14, 15, 16, 17, 18, 19)], Total_Bead_P1_K[, c(2, 14, 15, 16, 17, 18, 19)])
TotalPAR1_KO <- TotalPAR1_KO[ , -c(8)]

##Total Datasheet #PAR4-KO-----------------
Total_Input_P4_K <- read.csv("C:/Users/rahul/OneDrive/Desktop/HISAT_COMPARE/06-DE_lists/deseq2-Par4_ECko_Saline_Input_vs_Par4_ECko_APAP_Input_w_gene_names.tsv", sep = "")
Total_Bead_P4_K <- read.csv("C:/Users/rahul/OneDrive/Desktop/HISAT_COMPARE/06-DE_lists/deseq2-Par4_ECko_Saline_Bead_vs_Par4_ECko_APAP_Bead_w_gene_names.tsv", sep = "")

TotalPAR4_KO <- cbind(Total_Input_P4_K[, c(2, 14, 15, 16, 17, 18, 19)], Total_Bead_P4_K[, c(2, 14, 15, 16, 17, 18, 19)])
TotalPAR4_KO <- TotalPAR4_KO[ , -c(8)]

Total_K <- cbind(TotalPAR1_KO, TotalPAR4_KO)
Total_K <- Total_K[ , -c(14)]

##CombinedSheetsMeans----------------------

TotalSheet1 <- mutate(TotalPAR1_HET, Par1_EChet_Saline_Input = rowMeans(select(TotalPAR1_HET, X03.bam.Par1_EChet_Saline_Input_1.bam:X03.bam.Par1_EChet_Saline_Input_3.bam), na.rm = TRUE))
TotalSheet2 <- mutate(TotalPAR1_HET, Par1_EChet_APAP_Input = rowMeans(select(TotalPAR1_HET, X03.bam.Par1_EChet_APAP_Input_1.bam:X03.bam.Par1_EChet_APAP_Input_3.bam), na.rm = TRUE))
TotalSheet3 <- mutate(TotalPAR1_HET, Par1_EChet_Saline_Bead = rowMeans(select(TotalPAR1_HET, X03.bam.Par1_EChet_Saline_Bead_1.bam:X03.bam.Par1_EChet_Saline_Bead_3.bam), na.rm = TRUE))
TotalSheet4 <- mutate(TotalPAR1_HET, Par1_EChet_APAP_Bead = rowMeans(select(TotalPAR1_HET, X03.bam.Par1_EChet_APAP_Bead_1.bam:X03.bam.Par1_EChet_APAP_Bead_3.bam), na.rm = TRUE))

P1_H_Means <- cbind(TotalSheet1[, c(1,14)],TotalSheet2[, c(14)],TotalSheet3[, c(14)],TotalSheet4[, c(14)])

TotalSheet5 <- mutate(TotalPAR4_HET, Par4_EChet_Saline_Input = rowMeans(select(TotalPAR4_HET, X03.bam.Par4_EChet_Saline_Input_1.bam:X03.bam.Par4_EChet_Saline_Input_3.bam), na.rm = TRUE))
TotalSheet6 <- mutate(TotalPAR4_HET, Par4_EChet_APAP_Input = rowMeans(select(TotalPAR4_HET, X03.bam.Par4_EChet_APAP_Input_1.bam:X03.bam.Par4_EChet_APAP_Input_3.bam), na.rm = TRUE))
TotalSheet7 <- mutate(TotalPAR4_HET, Par4_EChet_Saline_Bead = rowMeans(select(TotalPAR4_HET, X03.bam.Par4_EChet_Saline_Bead_1.bam:X03.bam.Par4_EChet_Saline_Bead_3.bam), na.rm = TRUE))
TotalSheet8 <- mutate(TotalPAR4_HET, Par4_EChet_APAP_Bead = rowMeans(select(TotalPAR4_HET, X03.bam.Par4_EChet_APAP_Bead_1.bam:X03.bam.Par4_EChet_APAP_Bead_3.bam), na.rm = TRUE))

P4_H_Means <- cbind(TotalSheet5[, c(1,14)],TotalSheet6[, c(14)],TotalSheet7[, c(14)],TotalSheet8[, c(14)])
Het_Means <- cbind(P1_H_Means, P4_H_Means[, c(2,3,4,5)])


TotalSheet9 <- mutate(TotalPAR1_KO, Par1_ECko_Saline_Input = rowMeans(select(TotalPAR1_KO, X03.bam.Par1_ECko_Saline_Input_1.bam:X03.bam.Par1_ECko_Saline_Input_3.bam), na.rm = TRUE))
TotalSheet10 <- mutate(TotalPAR1_KO, Par1_ECko_APAP_Input = rowMeans(select(TotalPAR1_KO, X03.bam.Par1_ECko_APAP_Input_1.bam:X03.bam.Par1_ECko_APAP_Input_3.bam), na.rm = TRUE))
TotalSheet11 <- mutate(TotalPAR1_KO, Par1_ECko_Saline_Bead = rowMeans(select(TotalPAR1_KO, X03.bam.Par1_ECko_Saline_Bead_1.bam:X03.bam.Par1_ECko_Saline_Bead_3.bam), na.rm = TRUE))
TotalSheet12 <- mutate(TotalPAR1_KO, Par1_ECko_APAP_Bead = rowMeans(select(TotalPAR1_KO, X03.bam.Par1_ECko_APAP_Bead_1.bam:X03.bam.Par1_ECko_APAP_Bead_3.bam), na.rm = TRUE))

P1_K_Means <- cbind(TotalSheet9[, c(1,14)],TotalSheet10[, c(14)],TotalSheet11[, c(14)],TotalSheet12[, c(14)])

TotalSheet13 <- mutate(TotalPAR4_KO, Par4_ECko_Saline_Input = rowMeans(select(TotalPAR4_KO, X03.bam.Par4_ECko_Saline_Input_1.bam:X03.bam.Par4_ECko_Saline_Input_3.bam), na.rm = TRUE))
TotalSheet14 <- mutate(TotalPAR4_KO, Par4_ECko_APAP_Input = rowMeans(select(TotalPAR4_KO, X03.bam.Par4_ECko_APAP_Input_1.bam:X03.bam.Par4_ECko_APAP_Input_3.bam), na.rm = TRUE))
TotalSheet15 <- mutate(TotalPAR4_KO, Par4_ECko_Saline_Bead = rowMeans(select(TotalPAR4_KO, X03.bam.Par4_ECko_Saline_Bead_1.bam:X03.bam.Par4_ECko_Saline_Bead_3.bam), na.rm = TRUE))
TotalSheet16 <- mutate(TotalPAR4_KO, Par4_ECko_APAP_Bead = rowMeans(select(TotalPAR4_KO, X03.bam.Par4_ECko_APAP_Bead_1.bam:X03.bam.Par4_ECko_APAP_Bead_3.bam), na.rm = TRUE))

P4_K_Means <- cbind(TotalSheet9[, c(1,14)],TotalSheet10[, c(14)],TotalSheet11[, c(14)],TotalSheet12[, c(14)])

KO_Means <- cbind(P1_K_Means, P4_K_Means[, c(2,3,4,5)])

Means <- cbind(Het_Means, KO_Means[, c(2:9)])

title <- c("gene_name", "Par1_EChet_Saline_Input", "Par1_EChet_APAP_Input", "Par1_EChet_Saline_Bead", "Par1_EChet_APAP_Bead",
           "Par4_EChet_Saline_Input", "Par4_EChet_APAP_Input", "Par4_EChet_Saline_Bead", "Par4_EChet_APAP_Bead",
           "Par1_ECko_Saline_Input", "Par1_ECko_APAP_Input", "Par1_ECko_Saline_Bead", "Par1_ECko_APAP_Bead",
           "Par4_ECko_Saline_Input", "Par4_ECko_APAP_Input", "Par4_ECko_Saline_Bead", "Par4_ECko_APAP_Bead")

colnames(Means) <- title


##CombinedDataSheets Total----------------

TotalFinal <- cbind(Total, Total_K[, -c(1)])

ToTStr1 <- c("Par1 EChet Saline Input 1","Par1 EChet Saline Input 2","Par1 EChet Saline Input 3",
             "Par1 EChet APAP Input 1","Par1 EChet APAP Input 2","Par1 EChet APAP Input 3",
             "Par1 EChet Saline Bead 1","Par1 Saline EChet Bead 2","Par1 EChet Saline Bead 3",
             "Par1 EChet APAP Bead 1","Par1 EChet APAP Bead 2","Par1 EChet APAP Bead 3", 
             "Par4 EChet Saline Input 1","Par4 EChet Saline Input 2","Par4 EChet Saline Input 3",
             "Par4 EChet APAP Input 1","Par4 EChet APAP Input 2","Par4 EChet APAP Input 3",
             "Par4 EChet Saline Bead 1","Par4 EChet Saline Bead 2","Par4 EChet Saline Bead 3",
             "Par4 EChet APAP Bead 1","Par4 EChet APAP Bead 2","Par4 EChet APAP Bead 3")

ToTStr2 <- c("Par1 ECko Saline Input 1","Par1 ECko Saline Input 2","Par1 ECko Saline Input 3",
             "Par1 ECko APAP Input 1","Par1 ECko APAP Input 2","Par1 ECko APAP Input 3",
             "Par1 ECko Saline Bead 1","Par1 ECko Saline Bead 2","Par1 ECko Saline Bead 3",
             "Par1 ECko APAP Bead 1","Par1 ECko APAP Bead 2","Par1 ECko APAP Bead 3", 
             "Par4 ECko Saline Input 1","Par4 ECko Saline Input 2","Par4 ECko Saline Input 3",
             "Par4 ECko APAP Input 1","Par4 ECko APAP Input 2","Par4 APAP ECko Input 3",
             "Par4 ECko Saline Bead 1","Par4 ECko Saline Bead 2","Par4 ECko Saline Bead 3",
             "Par4 ECko APAP Bead 1","Par4 ECko APAP Bead 2","Par4 ECko APAP Bead 3")

colnames(TotalFinal) <- c(c("gene_name"), ToTStr1, ToTStr2)

#write.csv(x = TotalFinal, file = "FinalDataSheet.csv")
##------------------------------------------------------


##------------------------------------------------------
##Total Heatmap/PCA Plots #PAR1-HET-------------------
ToTStr3 <- c("Par1 EChet Saline Input 1","Par1 EChet Saline Input 2","Par1 EChet Saline Input 3",
            "Par1 EChet APAP Input 1","Par1 EChet APAP Input 2","Par1 EChet APAP Input 3",
            "Par1 EChet Saline Bead 1","Par1 EChet Saline Bead 2","Par1 EChet Saline Bead 3",
            "Par1 EChet APAP Bead 1","Par1 EChet APAP Bead 2","APAP Bead 3")
PrePro_heatprocess(TotalPAR1_HET, TotalPAR1_HET$gene_name, ToTStr3)
PrePro_PCAprocess(TotalPAR1_HET, TotalPAR1_HET$gene_name, ToTStr3, FALSE)

##Total Heatmap #PAR4-HET------------------------------
ToTStr4 <- c("Par4 EChet Saline Input 1","Par4 EChet Saline Input 2","Par4 EChet Saline Input 3",
            "Par4 EChet APAP Input 1","Par4 EChet APAP Input 2","Par4 EChet APAP Input 3",
            "Par4 EChet Saline Bead 1","Par4 EChet Saline Bead 2","Par4 EChet Saline Bead 3",
            "Par4 EChet APAP Bead 1","Par4 EChet APAP Bead 2","Par4 EChet APAP Bead 3")
PrePro_heatprocess(TotalPAR4_HET, TotalPAR4_HET$gene_name, ToTStr4)
PrePro_PCAprocess(TotalPAR4_HET, TotalPAR4_HET$gene_name, ToTStr4, FALSE)


PrePro_PCAprocess(Total, Total$gene_name, ToTStr1, FALSE)

##Total Heatmap #PAR1-KO------------------
ToTStr5 <- c("Par1 ECko Saline Input 1","Par1 ECko Saline Input 2","Par1 ECko Saline Input 3",
            "Par1 ECko APAP Input 1","Par1 ECko APAP Input 2","Par1 ECko APAP Input 3",
            "Par1 ECko Saline Bead 1","Par1 ECko Saline Bead 2","Par1 ECko Saline Bead 3",
            "Par1 ECko APAP Bead 1","Par1 ECko APAP Bead 2","Par1 ECko APAP Bead 3")
PrePro_heatprocess(TotalPAR1_KO, TotalPAR1_KO$gene_name, ToTStr5)
PrePro_PCAprocess(TotalPAR1_KO, TotalPAR1_KO$gene_name, ToTStr5, FALSE)

##Total Heatmap #PAR4-KO-----------------
ToTStr6 <- c("Par4 ECko Saline Input 1","Par4 ECko Saline Input 2","Par4 ECko Saline Input 3",
            "Par4 ECko APAP Input 1","Par4 ECko APAP Input 2","Par4 ECko APAP Input 3",
            "Par4 ECko Saline Bead 1","Par4 ECko Saline Bead 2","Par4 ECko Saline Bead 3",
            "Par4 ECko APAP Bead 1","Par4 ECko APAP Bead 2","Par4 ECko APAP Bead 3")
PrePro_heatprocess(TotalPAR4_KO, TotalPAR4_KO$gene_name, ToTStr6)
PrePro_PCAprocess(TotalPAR4_KO, TotalPAR4_KO$gene_name, ToTStr6, FALSE)


PrePro_PCAprocess(Total_K, Total_K$gene_name, ToTStr2, FALSE)


##Combined Heatmap and PCA Plot (Means)-----------------

title1 <- c("Par1_EChet_Saline_Input", "Par1_EChet_APAP_Input", "Par1_EChet_Saline_Bead", "Par1_EChet_APAP_Bead",
            "Par4_EChet_Saline_Input", "Par4_EChet_APAP_Input", "Par4_EChet_Saline_Bead", "Par4_EChet_APAP_Bead",
            "Par1_ECko_Saline_Input", "Par1_ECko_APAP_Input", "Par1_ECko_Saline_Bead", "Par1_ECko_APAP_Bead",
            "Par4_ECko_Saline_Input", "Par4_ECko_APAP_Input", "Par4_ECko_Saline_Bead", "Par4_ECko_APAP_Bead")


PrePro_heatprocess(Means, Means$gene_name, title1)
PrePro_PCAprocess(Means, Means$gene_name, title1, FALSE)

##Combined Heatmap and PCA Plot (Totals)-----------------
PrePro_heatprocess(TotalFinal, TotalFinal$gene_name, c(ToTStr1,ToTStr2))
PrePro_PCAprocess(TotalFinal, TotalFinal$gene_name, c(ToTStr1,ToTStr2), FALSE)


raw_TSheet <- TotalFinal
SIG_TSheet <- TotalFinal$gene_name
FR_K1 <- raw_TSheet[raw_TSheet$gene_name %in% SIG_TSheet, ]
j <- FR_K1$gene_name
FR_K1T <- as.data.frame(t(FR_K1[, -1]))
colnames(FR_K1T) <- j
FR_K1T$myfactor <- factor(row.names(FR_K1T))
colnames(FR_K1T) <- j
rownames(FR_K1T) <- c(ToTStr1,ToTStr2)
FR_K1T <- FR_K1T[, unlist(lapply(FR_K1T, is.numeric))]  # Remove non-numeric columns
FR_K1T <- FR_K1T[, !apply(FR_K1T == 0, 2, all)]         # Remove non zero columns
FR_K1T.pca <- prcomp(FR_K1T, center = TRUE, scale = TRUE)
summary(FR_K1T.pca)
str(FR_K1T.pca)
ggbiplot(FR_K1T.pca, choices = c(3,4), labels = rownames(FR_K1T), var.axes = FALSE)

prin_comp1 <- prcomp(FR_K1T, rank. = 3, scale = TRUE, center = TRUE)

ToTStr14 <- c("Par1 EChet Saline Input","Par1 EChet Saline Input","Par1 EChet Saline Input",
             "Par1 EChet APAP Input","Par1 EChet APAP Input","Par1 EChet APAP Input",
             "Par1 EChet Saline Bead","Par1 EChet Saline Bead","Par1 EChet Saline Bead",
             "Par1 EChet APAP Bead","Par1 EChet APAP Bead","Par1 EChet APAP Bead", 
             "Par4 EChet Saline Input","Par4 EChet Saline Input","Par4 EChet Saline Input",
             "Par4 EChet APAP Input","Par4 EChet APAP Input","Par4 EChet APAP Input",
             "Par4 EChet Saline Bead","Par4 EChet Saline Bead","Par4 EChet Saline Bead",
             "Par4 EChet APAP Bead","Par4 EChet APAP Bead","Par4 EChet APAP Bead")

ToTStr24 <- c("Par1 ECko Saline Input","Par1 ECko Saline Input","Par1 ECko Saline Input",
             "Par1 ECko APAP Input","Par1 ECko APAP Input","Par1 ECko APAP Input",
             "Par1 ECko Saline Bead","Par1 ECko Saline Bead","Par1 ECko Saline Bead",
             "Par1 ECko APAP Bead","Par1 ECko APAP Bead","Par1 ECko APAP Bead", 
             "Par4 ECko Saline Input","Par4 ECko Saline Input","Par4 ECko Saline Input",
             "Par4 ECko APAP Input","Par4 ECko APAP Input","Par4 ECko APAP Input",
             "Par4 ECko Saline Bead","Par4 ECko Saline Bead","Par4 ECko Saline Bead",
             "Par4 ECko APAP Bead","Par4 ECko APAP Bead","Par4 ECko APAP Bead")

components <- prin_comp1[["x"]]
components <- data.frame(components)
components$PC2 <- -components$PC2
components$PC3 <- -components$PC3
components = cbind(components, c(ToTStr14,ToTStr24))
colnames(components) <- c("PC1", "PC2", "PC3", "Sample")


tot_explained_variance_ratio <- summary(prin_comp1)[["importance"]]['Proportion of Variance',]
tot_explained_variance_ratio <- 100 * sum(tot_explained_variance_ratio)

tit = 'PCA Plot of PAR APAP Samples'

fig <- plot_ly(components, x = ~PC1, y = ~PC2, z = ~PC3, color = components$Sample, colors = c('#000000','#CD0BBC','#28E2E5','#61D04F',
                                                                                               '#018571','#5E3C99','#E66101','#A6611A',
                                                                                               '#CA0020','#003300','#0033FF','#ABA300',
                                                                                               '#00A9FF','#CCCFFF','#FF3300','#FFFF66'), size = 1000 ) %>%add_markers(size = 4)

fig <- fig %>%
  layout(
    title = tit,
    scene = list(bgcolor = "#FFFFFF")
  )

fig





PrePro_PCAprocess <- function(Total, GeneOI, Titles, AXES)
{
  Total <- TotalFinal
  GeneOI <- TotalFinal$gene_name
  Titles <- c(ToTStr1,ToTStr2)
  
  FR_K1 <- Total
  j <- FR_K1[FR_K1$gene_name %in% GeneOI,]
  FR_K1T <- as.data.frame(t(FR_K1[,-1]))
  colnames(FR_K1T) <- j
  FR_K1T$myfactor <- factor(row.names(FR_K1T))
  colnames(FR_K1T) <- j
  rownames(FR_K1T) <- Titles
  FR_K1T <- FR_K1T[ , unlist(lapply(FR_K1T, is.numeric))]# Remove non-numeric columns
  FR_K1T <- FR_K1T[, !apply(FR_K1T == 0, 2, all)] # Remove non zero columns
  FR_K1T.pca <-prcomp(FR_K1T, center = TRUE, scale = TRUE)
  summary(FR_K1T.pca)
  str(FR_K1T.pca)
  ggbiplot(FR_K1T.pca, labels=rownames(FR_K1T), var.axes = FALSE)
}


##------------------------------------------------------

##------------------------------------------------------
##Other Data Sheets-------------------------------------
##Cleaned APAP Bead fractions of PAR4 KO Bead (Removed Depleted Genes from APAP Liver)
View(dP4HAIvP4HAB_SIG)
dep_dP4HAIvP4HAB <- subset(dP4HAIvP4HAB_SIG, log2FoldChange<0)
dP4HABvP4KAB_clean <- dP4HABvP4KAB[ ! dP4HAIvP4HAB$gene_name %in% dep_dP4HAIvP4HAB$gene_name,]
dP4HABvP4KAB_clean$PAdj <- p.adjust(dP4HABvP4KAB_clean$PValue, method = "BH", n = nrow(dP4HABvP4KAB_clean))
EnhancedVolcano(dP4HABvP4KAB_clean, lab = dP4HABvP4KAB_clean$gene_name , x = 'log2FoldChange', y = 'PAdj', FCcutoff = 0.5, pCutoff = 0.05, title = "Par4 EChet Saline Bead vs Par4 EChet APAP Bead")

##Cleaned APAP Bead fractions of PAR4 KO Bead(Removed Saline Depleted Genes from Saline Liver)
View(dP4HSIvP4HSB_SIG)
dep_dP4HSIvP4HSB <- subset(dP4HSIvP4HSB_SIG, log2FoldChange<0)
dP4HABvP4KAB_clean <- dP4HABvP4KAB[ ! dP4HABvP4KAB$gene_name %in% dep_dP4HSIvP4HSB$gene_name,]
dP4HABvP4KAB_clean$PAdj <- p.adjust(dP4HABvP4KAB_clean$PValue, method = "BH", n = nrow(dP4HABvP4KAB_clean))
EnhancedVolcano(dP4HABvP4KAB_clean, lab = dP4HABvP4KAB_clean$gene_name , x = 'log2FoldChange', y = 'PAdj', FCcutoff = 0.5, pCutoff = 0.05, title = "Par4 EChet Saline Bead vs Par4 EChet APAP Bead")

##Cleaned APAP Bead fractions of PAR1 KO Bead(Removed Saline Depleted Genes from Saline Liver)
dep_dP1HSIvP1HSB <- subset(dP1HSIvP1HSB_SIG, log2FoldChange<0)
dP1HABvP1KAB_clean <- dP1HABvP1KAB[ ! dP1HABvP1KAB$gene_name %in% dep_dP1HSIvP1HSB$gene_name,]
dP1HABvP1KAB_clean$PAdj <- p.adjust(dP1HABvP1KAB_clean$PValue, method = "BH", n = nrow(dP1HABvP1KAB_clean))
dP1HABvP1KAB_clean_SIG <- subset(dP1HABvP1KAB_clean, PAdj <0.05)



EnhancedVolcano(dP1HABvP1KAB_clean, lab = dP1HABvP1KAB_clean$gene_name , x = 'log2FoldChange', y = 'PAdj', FCcutoff = 0.5, pCutoff = 0.05, title = "Par1 EChet Saline Bead vs Par1 EChet APAP Bead")
write.csv(x = dP1HABvP1KAB_clean, file = "cleaned PAr1 ecko bead apap.csv")


##Cleaned APAP Bead fractions of PAR1 EChet Bead Samples (Removed Saline Depleted Genes from Saline Liver)
dep_dP1HSIvP1HSB <- subset(dP1HSIvP1HSB_SIG, log2FoldChange<0)
dP1HSBvP1HAB_clean <- dP1HSBvP1HAB[ ! dP1HSBvP1HAB$gene_name %in% dep_dP1HSIvP1HSB$gene_name,]
dP1HSBvP1HAB_clean$PAdj <- p.adjust(dP1HSBvP1HAB_clean$PValue, method = "BH", n = nrow(dP1HABvP1KAB_clean))
dP1HSBvP1HAB_clean_SIG <- subset(dP1HSBvP1HAB_clean, PAdj <0.05)
write.csv(x = dP1HSBvP1HAB_clean_SIG, file = "cleaned PAr1 bead with apap_Sig.csv")
EnhancedVolcano(dP1HSBvP1HAB_clean, lab = dP1HSBvP1HAB_clean$gene_name , x = 'log2FoldChange', y = 'PAdj', FCcutoff = 0.5, pCutoff = 0.05, title = "Par1 EChet Saline Bead vs Par1 EChet APAP Bead")
write.csv(x = dP1HSIvP1HSB_SIG, file = "EC Enriched.csv")

##Cleaned APAP Bead fractions of PAR1 ECKO against PAR1 EChet Saline Bead Samples (Removed Saline Depleted Genes from Saline Liver)
dP1HSBvP1KAB <- D_subset(read.csv("C:/Users/rahul/OneDrive/Desktop/HISAT_COMPARE/06-DE_lists/deseq2-Par1_EChet_Saline_Bead_vs_Par1_ECko_APAP_Bead_w_gene_names.tsv", sep = ""), 10)
dP1HSBvP1KAB_clean <- dP1HSBvP1KAB[ ! dP1HSBvP1KAB$gene_name %in% dep_dP1HSIvP1HSB$gene_name,]
dP1HSBvP1KAB_clean$PAdj <- p.adjust(dP1HSBvP1KAB_clean$PValue, method = "BH", n = nrow(dP1HSBvP1KAB_clean))
dP1HSBvP1KAB_clean_SIG <- subset(dP1HSBvP1KAB_clean, PAdj <0.05)
write.csv(x = dP1HSBvP1KAB_clean, file = "dP1HSBvP1KAB_clean.csv")

##Cleaned Saline Bead fractions of PAR1 KO Bead(Removed Saline Depleted Genes from Saline Liver)
dP1HSBvP1KSB_clean <- dP1HSBvP1KSB[ ! dP1HSBvP1KSB$gene_name %in% dep_dP1HSIvP1HSB$gene_name,]
dP1HSBvP1KSB_clean$PAdj <- p.adjust(dP1HSBvP1KSB_clean$PValue, method = "BH", n = nrow(dP1HSBvP1KSB_clean))
dP1HSBvP1KSB_clean_SIG <- subset(dP1HSBvP1KSB_clean, PAdj <0.05)
EnhancedVolcano(dP1HSBvP1KSB_clean, lab = dP1HSBvP1KSB_clean$gene_name , x = 'log2FoldChange', y = 'PAdj', FCcutoff = 0.5, pCutoff = 0.05, title = "Par1 EChet Saline Bead vs Par1 ECko Saline Bead")

##------------------------------------------------------

##------------------------------------------------------
##Comparison Data of Altered Pathways in PAR1 ECKO Bead Samples------------------------------------------------------

cPath <- read.table("C:/Users/rahul/OneDrive/Desktop/HISAT_COMPARE/Comparisons_Pathway.txt", quote = "", sep = "\t", header=TRUE)
cPath <- cPath[!cPath$WT == "N/A", ]
cPath <- cPath[!cPath$KO == "N/A", ]
cPath <- cPath[!cPath$WT == 0, ]
cPath <- cPath[!cPath$KO == 0, ]

# Remove rows where column values fall between -0.5 and 0.5
#filtered_df <- cPath[!(cPath$WT >= -0.5 & cPath$WT <= 0.5), ]
#filtered_df <- cPath[!(cPath$KO >= -0.5 & cPath$KO <= 0.5), ]

filtered_df <- cPath

filtered_df$WT <- as.numeric(filtered_df$WT)
filtered_df$KO <- as.numeric(filtered_df$KO)
filtered_df$diff <- filtered_df$WT-filtered_df$KO

filtered_df$WT_CK <- ifelse(filtered_df$WT > 0, "POS", "NEG")
filtered_df$KO_CK <- ifelse(filtered_df$KO > 0, "POS", "NEG")

filtered_df$COM <- ifelse(filtered_df$WT_CK == filtered_df$KO_CK, "NO", "YES")
shifts <- sum(filtered_df$COM == "YES")
write.csv(x = filtered_df, file = "Comparisons Cleaned.csv")

##------------------------------------------------------


convert_numeric <- function(x) {
  if (!is.numeric(x)) {
    return(is.numeric(x))
  } else {
    return(x)
  }
}

# Call the function
TotalFinal[, -1] <- apply(TotalFinal[, -1], 2, convert_numeric)
TotalFinal$MEAN <- rowMeans(TotalFinal[, -1])


TotalFinal1 <- TotalFinal[TotalFinal$MEAN > 10, ]

##-----------------------------------------------------------------------

P1Data <- TotalFinal[, c(1, 8:13, 32:37)]

P1Data$EChet_Saline_Bead <- rowMeans(P1Data[, c(2,3,4)], na.rm = TRUE) 
P1Data$EChet_APAP_Bead <- rowMeans(P1Data[, c(5,6,7)], na.rm = TRUE) 
P1Data$ECko_Saline_Bead <- rowMeans(P1Data[, c(8,9,10)], na.rm = TRUE) 
P1Data$ECko_APAP_Bead <- rowMeans(P1Data[, c(11,12,13)], na.rm = TRUE) 

P1Data_Clean <- P1Data[, c(1, 14:17)]

PrePro_heatprocess(P1Data_Clean, c("F2r","Dll4","Arhgap15","Elmo1","Tuba1c","Tuba1b","Rims3","Arl4c","Pde7a","Ap2b1"), c("Par1 EChet Saline Bead","Par1 EChet APAP Bead","Par1 ECko Saline Bead", "Par1 ECko APAP Bead"))


