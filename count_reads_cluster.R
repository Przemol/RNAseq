require(GenomicAlignments, lib.loc = "/home/md1nbu/R/x86_64-pc-linux-gnu-library/3.3")
require(Rsamtools, lib.loc = "/home/md1nbu/R/x86_64-pc-linux-gnu-library/3.3")
require(SummarizedExperiment, lib.loc = "/home/md1nbu/R/x86_64-pc-linux-gnu-library/3.3")
require(TxDb.Hsapiens.UCSC.hg38.knownGene, lib.loc = "/home/md1nbu/R/x86_64-pc-linux-gnu-library/3.3")


gene_model <- exonsBy(keepStandardChromosomes(TxDb.Hsapiens.UCSC.hg38.knownGene), 'gene')
seqinfo(gene_model)

setwd('/mnt/fastdata/md1nbu/PRJNA369618/fastq')
bam_files <- BamFileList(dir(pattern = 'bam$'))
seqinfo(bam_files)

SE <- summarizeOverlaps(gene_model, bam_files)

str(assay(SE))

save(SE, file='SE.Rdata')