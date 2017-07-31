require(GenomicAlignments)
require(Rsamtools)
require(SummarizedExperiment)
require(TxDb.Hsapiens.UCSC.hg38.knownGene)


gene_model <- exonsBy(keepStandardChromosomes(TxDb.Hsapiens.UCSC.hg38.knownGene), 'gene')
seqinfo(gene_model)

setwd('/Users/przemol/BIGfiles/')
bam_files <- BamFileList(dir(pattern = 'bam$'))
seqinfo(bam_files)

SE <- summarizeOverlaps(gene_model, bam_files)

str(assay(SE))

