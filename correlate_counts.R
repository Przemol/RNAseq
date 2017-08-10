require(SummarizedExperiment)
require(DESeq2)
require(magrittr)


# Join Rdata files to singe SummarizedExperiment class
dir('data', full.names = TRUE, pattern = 'Rdata') %>% sapply(function(x) get(load(x)) ) -> SE_list
SE <- do.call(cbind, SE_list)
colnames(SE) <- gsub('(SRR[0-9]+)\\.bam', '\\1', colnames(SE))

# Extract count matrix
raw_counts <- assay(SE)

# Correalte
raw_counts <- raw_counts[!zeros, ] 
cor_mat <- cor(t(raw_counts))

cor_mat_no_diag <- cor_mat
diag(cor_mat_no_diag) <- 0

c1 <- cor_mat == 1
idx <- which(c1, arr.ind = TRUE)[1,]

# Asses significnace
library(psych)
cor_pvl <- corr.p(cor_mat , nrow(cor_mat), adjust="fdr", alpha=.05)


### translate ENTREZID to gene names
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
gene_smbl <- mapIds(org.Hs.eg.db, rownames(SE), "SYMBOL", "ENTREZID")
gene_names <- mapIds(org.Hs.eg.db, rownames(SE), "GENENAME", "ENTREZID")





# other way to acces counts
SE_list %>% sapply(assay) -> counts
colnames(counts) <- gsub('data/(SRR[0-9]+)\\.Rdata', '\\1', colnames(counts))







