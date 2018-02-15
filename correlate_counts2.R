require(SummarizedExperiment)
require(DESeq2)
require(magrittr)
library(org.Hs.eg.db)
library(edgeR)

# Join Rdata files to singe SummarizedExperiment class
setwd('/Users/nataliabulgakova/Documents/R/RNAseq')
dir('data', full.names = TRUE, pattern = 'Rdata') %>% sapply(function(x) get(load(x)) ) -> SE_list
SE <- do.call(cbind, SE_list)
colnames(SE) <- gsub('(SRR[0-9]+)\\.bam', '\\1', colnames(SE))

# Extract count matrix
raw_counts <- assay(SE)

# Correlate
raw_counts <- raw_counts[!zeros, ] 
raw_counts_cad <- raw_counts["999",]
cor_mat <- cor(t(raw_counts),raw_counts_cad, use='pairwise.complete.obs')


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
gene=c("CDH1")
unlist(mget(x=gene,envir=org.Hs.egALIAS2EG))




# other way to acces counts
SE_list %>% sapply(assay) -> counts
colnames(counts) <- gsub('data/(SRR[0-9]+)\\.Rdata', '\\1', colnames(counts))



###
load('cor_mat_rpkm_1.Rdata')

cor_mat2 <- cor_mat
printSig <- function(x=0.6) {
  library(tibble)
  library(dplyr)
  if(x>0){
    sig <- cor_mat2[cor_mat2>x,]
    tbl <- rev(sort(sig)) %>% as.data.frame() %>%  rownames_to_column('entrez') %>% tbl_df 
  } else {
    sig <- cor_mat2[cor_mat2<x,]
    tbl <- sort(sig) %>% as.data.frame() %>%  rownames_to_column('entrez') %>% tbl_df 
  }
  colnames(tbl) = c('entrez', 'cor')
  
  tbl$SYMBOL <- mapIds(org.Hs.eg.db, tbl$entrez, "SYMBOL", "ENTREZID")
  tbl$GENENAME <- mapIds(org.Hs.eg.db, tbl$entrez, "GENENAME", "ENTREZID")
  return(tbl)
}

a <- printSig(.5)
b <- printSig(-.5)

sig <- rbind(a,b)

cor_mat <- cor(t(assay(SE)[sig$entrez,]))

require(psych)
cor_pvl2 <- corr.p(cor_mat, nrow(assay(SE)), adjust="fdr", alpha=.05)


sort() %>% names %>%  
sort(cor_mat2[cor_mat2>0.6,]) %>% names %>%  mapIds(org.Hs.eg.db, ., "SYMBOL", "ENTREZID") %>% rev() %>% View()


# Correlate with STATs
raw_counts_cad <- raw_counts["999",]
cad <- t(t(raw_counts_cad))
raw_counts_stat <- raw_counts[c("6772", "6773", "6774", "6775", "6776", "6777", "6778"),]
statdata <- t(t(raw_counts_stat))

cor_matCS <- cor(t(raw_counts_stat),raw_counts_cad, use='pairwise.complete.obs')
cor_pvlCS <- corr.p(cor_matCS , nrow(cor_matCS), adjust="fdr", alpha=.05)

# Correlate with HP1 peaks
raw_counts_cad <- raw_counts["999",]
cad <- t(t(raw_counts_cad))
ortho <- read.csv(file="/Volumes/NGSDATA/sam/larvae/genes0-200orthologues2.csv", header=TRUE, sep=",")
ortho <- as.character(na.omit(ortho$Human.GeneID))
ortho <- unique(ortho)
ortho <- intersect(ortho, rownames(raw_counts))

raw_counts_ortho <- raw_counts[ortho,]
orthodata <- t(t(raw_counts_ortho))

cor_matortho <- cor(t(raw_counts_ortho),raw_counts_cad, use='pairwise.complete.obs')
cor_pvlortho <- corr.p(cor_matortho , nrow(cor_matortho), adjust="fdr", alpha=.05)

cor_matortho2 <- cor_matortho
library(tibble)
library(dplyr)
tbl2 <- as.data.frame(cbind(rownames(cor_matortho2),cor_pvlortho$r,cor_pvlortho$p))
colnames(tbl2) = c('entrez', 'cor','p-value')

tbl2$SYMBOL <- mapIds(org.Hs.eg.db, as.character(tbl2$entrez), "SYMBOL", "ENTREZID")
tbl2$GENENAME <- mapIds(org.Hs.eg.db, as.character(tbl2$entrez), "GENENAME", "ENTREZID")
tbl2 <- tbl2[tbl2$'p-value'<0.01,]


sort() %>% names %>%  
  sort(cor_pvlortho) %>% names %>%  mapIds(org.Hs.eg.db, ., "SYMBOL", "ENTREZID") %>% rev() %>% View()


# Remove low expressed genes
myCPM <- cpm(raw_counts)
thresh <- myCPM > 0.1
keep <- rowSums(thresh) >= 2
# plot(myCPM[,1],raw_counts[,1],ylim=c(0,50),xlim=c(0,3))
# abline(h=10)
counts.keep <- raw_counts[keep,]

# Convert counts to DGEList object
countsdata <- DGEList(counts.keep)
countsdata <- calcNormFactors(countsdata)
norm.counts <- counts.keep %*% diag(countsdata$samples$norm.factors)
norm_counts_stat <- t(norm.counts[c("6772", "6773", "6774", "6775", "6776", "6777", "6778"),])
