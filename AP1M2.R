if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")

library(TCGAbiolinks)
library(dplyr)
library(DT)
library(SummarizedExperiment)
require(psych)

query <- GDCquery(project = "TCGA-UCEC",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "STAR - Counts")


GDCdownload(query)
data <- GDCprepare(query)

raw_counts <- assay(data)
tab <- rowData(data)
map <- tab$gene_name
names(map) <- tab$gene_id
rownames(raw_counts) <- map[rownames(raw_counts)]

raw_counts_GOI <- raw_counts["CTNND1",]
raw_counts_Act <- raw_counts["ACTB",]

cor_mat <- cor(t(raw_counts),raw_counts_GOI, use='pairwise.complete.obs')
cor_mat_p <- corr.p(t(raw_counts),raw_counts_GOI , adjust="fdr", alpha=.05)

Positive <- sort(cor_mat[!is.na(cor_mat[,1]) & cor_mat>0.4,], decreasing = TRUE)
Negative <- sort(cor_mat[!is.na(cor_mat[,1]) & cor_mat< - 0.4,], decreasing = FALSE)
