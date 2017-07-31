require(GenomicAlignments)
require(Rsamtools)
require(SummarizedExperiment)
require(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)


gene_model <- exonsBy(keepStandardChromosomes(TxDb.Hsapiens.UCSC.hg38.knownGene), 'gene')
seqinfo(gene_model)

setwd('/Users/przemol/BIGfiles/')
bam_files <- BamFileList(dir(pattern = 'bam$'))
seqinfo(bam_files)

SE <- summarizeOverlaps(gene_model, bam_files)
str(assay(SE))


library(TxDb.Hsapiens.UCSC.hg38.knownGene)
gene_smbl <- mapIds(org.Hs.eg.db, rownames(SE), "SYMBOL", "ENTREZID")
gene_names <- mapIds(org.Hs.eg.db, rownames(SE), "GENENAME", "ENTREZID")


