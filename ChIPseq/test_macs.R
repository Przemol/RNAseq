library(rtracklayer)
library(Rsamtools)
library(GenomicRanges)
library(IRanges)

source("https://bioconductor.org/biocLite.R")
biocLite("org.Dm.eg.db")
biocLite("chipseq")
biocLite("BayesPeak")
library(org.Dm.eg.db)
setwd('/Volumes/NGSDATA/sam')

# Load and tile genome

biocLite("BSgenome.Dmelanogaster.UCSC.dm6")
library(BSgenome.Dmelanogaster.UCSC.dm6)
genome <- BSgenome.Dmelanogaster.UCSC.dm6
si <- seqinfo(genome)
library(biomaRt)
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "dmelanogaster_gene_ensembl",
                host="www.ensembl.org")

ensembl <- useEnsembl(biomart="ensembl", dataset="dmelanogaster_gene_ensembl")

ds <- useDataset('dmelanogaster_gene_ensembl', mart=mart)
egs <- getBM(attributes = c('ensembl_gene_id',
                            'chromosome_name','start_position',
                            'end_position','strand', 'external_gene_name', 'go_id'),mart=ds)
binsize <- 200
bins <- tileGenome(si, tilewidth=binsize, cut.last.tile.in.chrom=TRUE)

egs$TSS <- ifelse( egs$strand == "1", egs$start_position, egs$end_position )
head(egs)
promoter_size <- 200
promoter_regions <-
  GRanges(seqnames = Rle( paste0('chr', egs$chromosome_name) ),
          ranges = IRanges( start = egs$TSS - promoter_size,
                            end = egs$TSS + promoter_size ),
          strand = Rle( rep("*", nrow(egs)) ),
          gene = egs$ensembl_gene_id,
          	 name = egs$external_gene_name, go = egs$go_id)


# Load reads
input1 <- import.bed('3_Input_2228_filtered.bed')
input2 <- import.bed('4_Input_2229_filtered.bed')

rep1 <- import.bed('3_HP1a_2265_filtered.bed')
rep2 <- import.bed('4_HP1a_2266_filtered.bed')

# Prepare samples
library(chipseq)
prepareChIPseq <- function(reads){
  frag.len = median( estimate.mean.fraglen(reads), na.rm=TRUE )
  cat( paste0( 'Median fragment size for this library is ', round(frag.len)))
  reads.extended = resize(reads, width = frag.len)
  return( trim(reads.extended) )
}

input1 <- prepareChIPseq( input1 )
input2 <- prepareChIPseq( input2 )
rep1 <- prepareChIPseq( rep1 )
rep2 <- prepareChIPseq( rep2 )

# Binning
BinChIPseq = function( reads, bins ){
  mcols(bins)$score = countOverlaps( bins, reads )
  6
  return( bins )
}

input1.200bins = BinChIPseq( input1, bins )
input2.200bins = BinChIPseq( input2, bins )
rep1.200bins = BinChIPseq( rep1, bins )
rep2.200bins = BinChIPseq( rep2, bins )

# Identification of peaks
library(BayesPeak)
raw.output3 <- bayespeak('3_HP1a_2265_filtered.bed', '3_Input_2228_filtered.bed')
output3 <- summarize.peaks(raw.output3)
export.bed(output3, name = 'larva3')
raw.output4 <- bayespeak('4_HP1a_2266_filtered_ucsc.bed', '4_Input_2229_filtered_ucsc.bed')
export.bed(raw.output4, name = 'larva4')

peaks.rep1 <- import.bed('Rep1_summits.bed')
peaks.rep2 <- import.bed('Rep2_summits.bed')

# Overlap promoter with peaks

regions1 <- unique(subsetByOverlaps(promoter_regions,peaks.rep1))
regions2 <- unique(subsetByOverlaps(promoter_regions, peaks.rep2))

library(GenomicRanges)
library(data.table)

reg1 <- as.data.table(regions1)
reg2 <- as.data.table(regions2)
library(plyr)
genes <- join_all(list(reg1,reg2), by = 'name', type = 'inner')
genes <- genes[,6:8]

write.csv( genes, file = 'result.csv')
