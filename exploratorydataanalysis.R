##########################################################################################################
#---------------------------------------------------------------------------------------------------------
##########################################################################################################
#Exploratory Data Analysis

#set directory
setwd(dir = "C:/Users/Orian/OneDrive/Documents/Master")

#install packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install('EnhancedVolcano')
BiocManager::install("DESeq2")
BiocManager::install("biomaRt")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Mm.eg.db")


#load library
library(DESeq2)
library(biomaRt)
library(org.Mm.eg.db)
library(EnhancedVolcano)
library(clusterProfiler)


#import the count table and the column table (with information about the samples)
count_table <- read.table("final_count_table.txt", header = T, row.names = 1, sep = "\t")
col_table <- read.table("data_table.txt", header = T, row.names = 1, sep = "\t")

#create the DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = count_table, colData = col_table, design = ~ Group)

#perform the DE
dds <- DESeq(dds)

#PCA to vizualise the data
vsd <- vst(dds, blind = TRUE) #to remove the variance on the mean 
plotPCA(vsd, intgroup=c("Group"))

##########################################################################################################
#---------------------------------------------------------------------------------------------------------
##########################################################################################################
#Differential Expression Analysis

###############################
# Blood Case vs Blood Control #
###############################

results_blood <- results(dds, contrast=c("Group","Blood_WT_Case","Blood_WT_Control"))
DE_blood <- results_blood[which(results_blood$padj < 0.05),]
summary(DE_blood)


#add genes name from ensembl, need to specify the version of the reference used previously (in our case, version 110)
ensembl <- useEnsembl(biomart = 'genes', dataset = 'mmusculus_gene_ensembl', version = 110)
ensembl_gene_ids_b <- rownames(DE_blood)
gene_info_b <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                   filters = "ensembl_gene_id",
                   values = ensembl_gene_ids_b,
                   mart = ensembl)

DE_blood_with_gene_names <- cbind(gene_info_b, DE_blood)

#sort DE genes by pvalue to get the ten most DE genes
sort_blood <- DE_blood_with_gene_names[order(DE_blood_with_gene_names$padj),]
sort_blood[1:10,]


#volcano plot to visualize DE genes
EnhancedVolcano(DE_blood_with_gene_names,
                lab = DE_blood_with_gene_names$external_gene_name,
                x = 'log2FoldChange',
                y = 'pvalue',
                col = c("pink3", "olivedrab", "cyan3", "coral1"),
                title = "DE genes - Blood Case vs Blood Control")

#explore expression of genes from the article
Rsad2 <-DE_blood["ENSMUSG00000020641",]
Rsad2$log2FoldChange > 0 #False --> underexpression
Rsad2 <- plotCounts(dds, "ENSMUSG00000020641", intgroup = c("Group"), returnData = TRUE)
boxplot(count ~ Group , data=Rsad2, main = "Expression of Rsad2")

oas1a <-DE_blood["ENSMUSG00000052776",]
oas1a$log2FoldChange > 0 #true -- >overexpression
oas1a <- plotCounts(dds, "ENSMUSG00000052776", intgroup = c("Group"), returnData = TRUE)
boxplot(count ~ Group , data=oas1a, main = "Expression of oas1a")

#############################
# Lung Case vs Lung Control #
#############################

results_lung <- results(dds, contrast=c("Group","Lung_WT_Case","Lung_WT_Control"))
DE_lung <- results_lung[which(results_lung$padj < 0.05),]
summary(DE_lung)

#add genes name
ensembl_gene_ids_l <- rownames(DE_lung)
gene_info_l <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                   filters = "ensembl_gene_id",
                   values = ensembl_gene_ids_l,
                   mart = ensembl)

DE_lung_with_gene_names <- cbind(gene_info_l, DE_lung)

sort_lung <- DE_lung_with_gene_names[order(DE_lung_with_gene_names$padj),]
sort_lung[1:10,]



#volcano plot
EnhancedVolcano(DE_lung_with_gene_names,
                lab = DE_lung_with_gene_names$external_gene_name,
                x = 'log2FoldChange',
                y = 'pvalue',
                col = c("pink3", "olivedrab", "cyan3", "coral1"),
                title = "DE genes -  Lung Case vs Lung Control")

#explore expression of genes from the article
irf7 <-DE_lung["ENSMUSG00000025498",]
irf7$log2FoldChange > 0 #True --> overexpression
irf7 <- plotCounts(dds, "ENSMUSG00000025498", intgroup = c("Group"), returnData = TRUE)
boxplot(count ~ Group , data=irf7, main = "Expression of irf7")

fcgr1 <-DE_lung["ENSMUSG00000015947",]
fcgr1$log2FoldChange > 0 #True -- >overexpression
fcgr1 <- plotCounts(dds, "ENSMUSG00000015947", intgroup = c("Group"), returnData = TRUE)
boxplot(count ~ Group , data=fcgr1, main = "Expression of fcgr1")

##########################################################################################################
#---------------------------------------------------------------------------------------------------------
##########################################################################################################
#Overexpression Analysis 

###############################
# Blood Case vs Blood Control #
###############################

go_blood <- enrichGO(gene = row.names(DE_blood_with_gene_names), universe = names(dds), 
                     OrgDb = org.Mm.eg.db, ont= "BP", keyType = "ENSEMBL")

head(go_blood)

#barplot: count by go terms, sorted by p-value
barplot(go_blood, showCategory = 15) + ggtitle("GO terms - Blood Case vs Blood Control")

#web plot: show relationship between go termes
goplot(go_blood, showCategory = 10) + ggtitle("GO terms - Blood Case vs Blood Control")

#dot plot: gene ratio
dotplot(go_blood) + ggtitle("GO terms - Blood Case vs Blood Control")

#############################
# Lung Case vs Lung Control #
#############################

go_lung <- enrichGO(gene = row.names(DE_lung_with_gene_names), universe = names(dds), 
                     OrgDb = org.Mm.eg.db, ont= "BP", keyType = "ENSEMBL")

head(go_lung)

#barplot: count by go terms, sorted by p-value
barplot(go_lung, showCategory = 15) + ggtitle("GO terms - Lung Case vs Lung Control")

#web plot: show relationship between go termes
goplot(go_lung, showCategory = 10) + ggtitle("GO terms - Lung Case vs Lung Control")

#dotplot: gene ratio
dotplot(go_lung) + ggtitle("GO terms - Lung Case vs Lung Control")

