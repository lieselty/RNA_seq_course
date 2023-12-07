#RNA-seq course, 09_exploratory
setwd(dir = "C:/Users/Orian/OneDrive/Documents/Master")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
library(DESeq2)
install.packages("pheatmap")
library(pheatmap)

count_table <- read.table("final_count_table.txt", header = T, row.names = 1, sep = "\t")
col_table <- read.table("data_table.txt", header = T, row.names = 1, sep = "\t")

dds <- DESeqDataSetFromMatrix(countData = count_table, colData = col_table, design = ~ Group)
dds <- DESeq(dds)

vsd <- vst(dds, blind = TRUE)
rld <- rlog(dds, blind = TRUE)

plotPCA(vsd, intgroup=c("Group"))

#differential expression analysis
#blood
results_blood <- results(dds, contrast=c("Group","Blood_WT_Case","Blood_WT_Control"), 
                         independentFiltering=FALSE, cooksCutoff=FALSE) #??
results_blood <- results_blood[which(results_blood$padj < 0.5),]
summary(results_blood) #is this enough?

Rsad2 <-results_blood["ENSMUSG00000020641",]
Rsad2$log2FoldChange > 0 #False --> underexpression
oas1a <-results_blood["ENSMUSG00000052776",]
oas1a$log2FoldChange > 0 #true -- >overexpression

##############################
#-----------------------------
##############################
#lung
results_lung <- results(dds, contrast=c("Group","Lung_WT_Case","Lung_WT_Control"),
                        independentFiltering=FALSE, cooksCutoff=FALSE)
results_lung <- results_lung[which(results_lung$padj < 0.5),]
summary(results_lung)

irf7 <-results_lung["ENSMUSG00000025498",]
irf7$log2FoldChange > 0 #True --> overexpression
fcgr1 <-results_lung["ENSMUSG00000015947",]
fcgr1$log2FoldChange > 0 #True -- >overexpression


