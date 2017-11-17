library(rtracklayer)
library(Rsamtools)
library(GenomicRanges)
library(IRanges)
library(data.table)
library(GenomicRanges)
source("https://bioconductor.org/biocLite.R")
biocLite("org.Dm.eg.db")
biocLite("chipseq")
biocLite("BayesPeak")
setwd('/Volumes/NGSDATA/sam/larvae')

# Load and tile genome

biocLite("BSgenome.Dmelanogaster.UCSC.dm6")
library(BSgenome.Dmelanogaster.UCSC.dm6)
genome <- BSgenome.Dmelanogaster.UCSC.dm6

library(biomaRt)
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "dmelanogaster_gene_ensembl",
                host="www.ensembl.org")
ds <- useDataset('dmelanogaster_gene_ensembl', mart=mart)
egs <- getBM(attributes = c('ensembl_gene_id',
                            'chromosome_name','start_position',
                            'end_position','strand', 'external_gene_name', 'go_id'),mart=ds)
head(egs)
egs$TSS <- ifelse( egs$strand == "1", egs$start_position, egs$end_position )
egs$upstream <- ifelse( egs$strand == "1", egs$TSS-500, egs$TSS+500)
promoter_size <- 200
promoter_regions <-
  GRanges(seqnames = Rle( paste0('chr', egs$chromosome_name) ),
          ranges = IRanges( start = egs$TSS - promoter_size,
                            end = egs$TSS + promoter_size ),
          strand = Rle( rep("*", nrow(egs)) ),
          gene = egs$ensembl_gene_id,
          name = egs$external_gene_name, go = egs$go_id)

promoter_regions2 <-
  GRanges(seqnames = Rle( paste0('chr', egs$chromosome_name) ),
          ranges = IRanges( start = ifelse(egs$upstream<egs$TSS,egs$upstream,egs$TSS),
                            end = ifelse(egs$upstream<egs$TSS,egs$TSS, egs$upstream) ),
          strand = Rle( rep("*", nrow(egs)) ),
          gene = egs$ensembl_gene_id,
          name = egs$external_gene_name, go = egs$go_id)

# Identification of peaks
library(BayesPeak)
raw.output3 <- bayespeak('3_HP1a_2265_filtered.bed', '3_Input_2228_filtered.bed')
output3 <- summarize.peaks(raw.output3, threshold = 0.8)
outdir3 <- file.path('/Volumes/NGSDATA/sam/larvae', "larva3.bed")
export.bed(output3,outdir3)

raw.output4 <- bayespeak('4_HP1a_2266_filtered.bed', '4_Input_2229_filtered.bed')
output4 <- summarize.peaks(raw.output4, threshold = 0.8)
outdir4 <- file.path('/Volumes/NGSDATA/sam/larvae', "larva4.bed")
export.bed(output4,outdir4)

enriched.regions <- Reduce(subsetByOverlaps, list(output3, output4))

enriched_regions <-
  GRanges(seqnames = Rle( enriched.regions$space ),
          ranges = IRanges( enriched.regions$ranges ),
          seqinfo=seqinfo(genome))

genes <- unique(subsetByOverlaps(promoter_regions,enriched_regions))
genes_df <- as(genes, "data.frame")
write.table(x = genes_df, 'genes200-200.csv', sep=",", col.names=TRUE)

genes2  <- unique(subsetByOverlaps(promoter_regions2,enriched_regions))
genes2_df <- as(genes2, "data.frame")
write.table(x = genes_df, 'genes0-500.csv', sep=",", col.names=TRUE)
