require(SummarizedExperiment, lib.loc = "/home/md1nbu/R/x86_64-pc-linux-gnu-library/3.3")
require(DESeq2, lib.loc = "/home/md1nbu/R/x86_64-pc-linux-gnu-library/3.3")
require(magrittr, lib.loc = "/home/md1nbu/R/x86_64-pc-linux-gnu-library/3.3")
library(psych, lib.loc = "/home/md1nbu/R/x86_64-pc-linux-gnu-library/3.3")

# Join Rdata files to singe SummarizedExperiment class
dir('/data/md1nbu/RNAseq/Data', full.names = TRUE, pattern = 'Rdata') %>% sapply(function(x) get(load(x)) ) -> SE_list
SE <- do.call(cbind, SE_list)
colnames(SE) <- gsub('(SRR[0-9]+)\\.bam', '\\1', colnames(SE))

# Extract count matrix
raw_counts <- assay(SE)
dds <- DESeqDataSet(SE, ~1)
rpkm <- fpkm(dds, robust = TRUE)
base_mean <- rowMeans(rpkm)
rpkm_flt <- rpkm[base_mean>1,]
save(rpkm_flt, file='cor_mat_rpkm_1.Rdata')

# Correalte
cor_mat <- cor(t(rpkm_flt))
save(cor_mat, file='cor_mat.Rdata')


# Asses significnace
cor_pvl <- corr.p(cor_mat , nrow(cor_mat), adjust="fdr", alpha=.05)
save(cor_pvl, file='cor_pvl.Rdata')





