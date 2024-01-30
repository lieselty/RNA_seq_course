#set directory
setwd(dir = "C:/Users/Orian/OneDrive/Documents/Master")

#install packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("biomaRt")
BiocManager::install("org.Mm.eg.db")
BiocManager::install('EnhancedVolcano')
BiocManager::install("clusterProfiler")

#load library
library(DESeq2)
library(biomaRt)
library(org.Mm.eg.db)
library(EnhancedVolcano)
library(clusterProfiler)


#############################
# Exploratory Data Analysis #------------------------------------------------------------------------------------------------------------
#############################

#import the count table and the column table (with information about the samples)
count_table <- read.table("final_count_table.txt", header = T, row.names = 1, sep = "\t")
col_table <- read.table("data_table.txt", header = T, row.names = 1, sep = "\t")


#create the DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = count_table, colData = col_table, design = ~ Group)

#perform the DE
dds <- DESeq(dds)


#PCA to vizualise the data
vsd <- vst(dds, blind = TRUE) #to remove the variance on the mean 
plotPCA(vsd, intgroup = c("Group")) + ggtitle("PCA") + theme(plot.title = element_text(hjust = 0.5))


####################################
# Differential Expression Analysis #-----------------------------------------------------------------------------------------------------
####################################

# # # # # # # # # # # # # # # #
# Blood Case vs Blood Control #
# # # # # # # # # # # # # # # #

results_blood <- results(dds, contrast=c("Group","Blood_WT_Case","Blood_WT_Control"))
#filter the results by padj < 0.05
DE_blood <- results_blood[which(results_blood$padj < 0.05),]
#To see number of DE genes:
summary(DE_blood)


#add genes name from ensembl
ensembl <- useEnsembl(biomart = 'genes', dataset = 'mmusculus_gene_ensembl', version = 110)#specify version
ensembl_gene_ids_b <- rownames(DE_blood)
gene_info_b <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                   filters = "ensembl_gene_id",
                   values = ensembl_gene_ids_b,
                   mart = ensembl)

DE_blood_with_gene_names <- cbind(gene_info_b, DE_blood)


#volcano plot to visualize DE genes
blood_volcano <- EnhancedVolcano(DE_blood_with_gene_names,
                lab = DE_blood_with_gene_names$external_gene_name,
                x = 'log2FoldChange',
                y = 'pvalue',
                col = c("pink3", "olivedrab", "cyan3", "coral1"),
                title = "DE genes - Blood Case vs Blood Control")


#explore expression of gene from the article
oas1a <-DE_blood["ENSMUSG00000052776",]
oas1a$log2FoldChange > 0 #true -- >overexpression in blood
oas1a <- plotCounts(dds, "ENSMUSG00000052776", intgroup = c("Group"), returnData = TRUE)
boxplot(count ~ Group , data=oas1a, main = "Expression of Oas1a")


# # # # # # # # # # # # # # #
# Lung Case vs Lung Control #
# # # # # # # # # # # # # # #

results_lung <- results(dds, contrast=c("Group","Lung_WT_Case","Lung_WT_Control"))
#filter the results by padj < 0.05
DE_lung <- results_lung[which(results_lung$padj < 0.05),]
#To see the number of DE genes:
summary(DE_lung)


#add genes name
ensembl_gene_ids_l <- rownames(DE_lung)
gene_info_l <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                   filters = "ensembl_gene_id",
                   values = ensembl_gene_ids_l,
                   mart = ensembl)

DE_lung_with_gene_names <- cbind(gene_info_l, DE_lung)


#volcano plot
lung_volcano <- EnhancedVolcano(DE_lung_with_gene_names,
                lab = DE_lung_with_gene_names$external_gene_name,
                x = 'log2FoldChange',
                y = 'pvalue',
                col = c("pink3", "olivedrab", "cyan3", "coral1"),
                title = "DE genes -  Lung Case vs Lung Control")


#explore expression of gene from the article
fcgr1 <-DE_lung["ENSMUSG00000015947",]
fcgr1$log2FoldChange > 0 #True -- >overexpression in lung
fcgr1 <- plotCounts(dds, "ENSMUSG00000015947", intgroup = c("Group"), returnData = TRUE)
boxplot(count ~ Group , data=fcgr1, main = "Expression of Fcgr1")



###########################
# Overexpression Analysis #--------------------------------------------------------------------------------------------------------------
###########################

# # # # # # # # # # # # # # # #
# Blood Case vs Blood Control #
# # # # # # # # # # # # # # # #

go_blood <- enrichGO(gene = row.names(DE_blood_with_gene_names), universe = names(dds), 
                     OrgDb = org.Mm.eg.db, ont= "BP", keyType = "ENSEMBL")


#dot plot: gene ratio
blood_dot <- dotplot(go_blood) + ggtitle("Blood (Case vs Control): GeneRatio")

#other plots not used in the paper
#barplot: count by go terms, sorted by p-value
blood_bar <- barplot(go_blood, showCategory = 10) + ggtitle("Blood (Case vs Control): P-value")

#web plot: show relationship between go terms
goplot(go_blood, showCategory = 10) + ggtitle("GO terms - Blood Case vs Blood Control")



# # # # # # # # # # # # # # #
# Lung Case vs Lung Control #
# # # # # # # # # # # # # # #

go_lung <- enrichGO(gene = row.names(DE_lung_with_gene_names), universe = names(dds), 
                     OrgDb = org.Mm.eg.db, ont= "BP", keyType = "ENSEMBL")


#dotplot: gene ratio
lung_dot <- dotplot(go_lung) + ggtitle("Lung (Case vs Control): GeneRatio")

#other plots not used in the paper
#barplot: count by go terms, sorted by p-value
lung_bar <- barplot(go_lung, showCategory = 10) + ggtitle("Lung (Case vs Control): P-value")

#web plot: show relationship between go terms
goplot(go_lung, showCategory = 10) + ggtitle("GO terms - Lung Case vs Lung Control")