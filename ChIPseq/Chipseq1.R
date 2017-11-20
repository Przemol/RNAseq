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
egs$upstream2 <- ifelse( egs$strand == "1", egs$TSS-200, egs$TSS+200)
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

promoter_regions3 <-
  GRanges(seqnames = Rle( paste0('chr', egs$chromosome_name) ),
          ranges = IRanges( start = ifelse(egs$upstream2<egs$TSS,egs$upstream2,egs$TSS),
                            end = ifelse(egs$upstream2<egs$TSS,egs$TSS, egs$upstream2) ),
          strand = Rle( rep("*", nrow(egs)) ),
          gene = egs$ensembl_gene_id,
          name = egs$external_gene_name, go = egs$go_id)

# Identification of peaks
library(BayesPeak)
raw.output3 <- bayespeak('3_HP1a_2265_filtered.bed', '3_Input_2228_filtered.bed')
output3 <- summarize.peaks(raw.output3, threshold = 0.9)
outdir3 <- file.path('/Volumes/NGSDATA/sam/larvae', "larva3.bed")
export.bed(output3,outdir3)

raw.output4 <- bayespeak('4_HP1a_2266_filtered.bed', '4_Input_2229_filtered.bed')
output4 <- summarize.peaks(raw.output4, threshold = 0.9)
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

genes3  <- unique(subsetByOverlaps(promoter_regions3,enriched_regions))
genes3_df <- as(genes3, "data.frame")
write.table(x = genes_df, 'genes0-200.csv', sep=",", col.names=TRUE)

# Gene selection by function and plotting
go.table <- read.csv('go_table.csv')

selected_genes <- go.table
selected_genes <- egs[egs$ensembl_gene_id %in% selected_genes$ID,] <- egs[egs$ensembl_gene_id %in% selected_genes$ID,]
selected_genes <- selected_genes[!duplicated(selected_genes$ensembl_gene_id),]


selected_genes$TSS <- ifelse( selected_genes$strand == "1", selected_genes$start_position, selected_genes$end_position )

tiles = sapply( 1:nrow(selected_genes), function(i)
  if( selected_genes$strand[i] == "1" )
    selected_genes$TSS[i] + seq( -1000, 900, length.out=20 )
  else
    selected_genes$TSS[i] + seq( 900, -1000, length.out=20 ) )

tiles = GRanges(tilename = paste( rep( selected_genes$ensembl_gene_id, each=20), 1:20, sep="_" ),
                seqnames = Rle( rep(paste0('chr', selected_genes$chromosome_name), each=20) ),
                ranges = IRanges(start = as.vector(tiles),
                                 width = 100),
                strand = Rle(rep("*", length(as.vector(tiles)))),
                seqinfo=seqinfo(genome))

# Load reads
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

rep1 <- prepareChIPseq( rep1 )
rep2 <- prepareChIPseq( rep2 )

# Count overlaps
HP1.p <- countOverlaps( tiles, rep1) +
  countOverlaps( tiles, rep2 )

dimnames1 <- list(selected_genes$ensembl_gene_id, 1:20)
HP1.p.matrix <- matrix( HP1.p, nrow=nrow(selected_genes),
                           ncol=20, byrow=TRUE, dimnames = dimnames1)

# Making graph

colors <- colorRampPalette(c('white','red','gray','black'))(100)
layout(mat=matrix(c(1,2,0,3), 2, 2),
       widths=c(2,2,2),
       heights=c(0.5,5,0.5,5), TRUE)
par(mar=c(2.5,1,1.5,1))
image(seq(0, max(HP1.p.matrix), length.out=100), 1,
      matrix(seq(0, max(HP1.p.matrix), length.out=100),100,1),
      col = colors,
      xlab='Distance from TSS', ylab='',
      main='Number of reads', yaxt='n',
      lwd=3, axes=TRUE)
box(col='black', lwd=2)
image(x=seq(-1000, 1000, length.out=20),
      y=1:nrow(HP1.p.matrix),
      z=t(HP1.p.matrix[order(rowSums(HP1.p.matrix)),]),
      col=colors,
      xlab='Distance from TSS (bp)',
      ylab='Promoters', lwd=2)
box(col='black', lwd=2)
abline(v=0, lwd=1, col='gray')

plot(x=seq(-1000, 1000, length.out=20),
     y=colMeans(HP1.p.matrix),
     ty='b', pch=19,
     col='red4',lwd=2,
     ylab='Mean tag count',
     xlab='Distance from TSS (bp)')
abline(h=seq(1,100,by=5),
       v=seq(-1000, 1000, length.out=20),
       lwd=0.25, col='gray')
box(col='black', lwd=2)

# Pre vs after TSS enrichment
HPafterTSS <- HP1.p.matrix[rowSums(HP1.p.matrix[,1:10])<rowSums(HP1.p.matrix[,11:20]),]
genesafterTSS <- rownames(HPafterTTS)
write.csv(genesafterTSS, file = 'genes_with_HP_after_TSS.csv')
HPpreTSS <- HP1.p.matrix[rowSums(HP1.p.matrix[,1:10])>rowSums(HP1.p.matrix[,11:20]),]
genespreTSS <- rownames(HPpreTTS)
write.csv(genespreTSS, file = 'genes_with_HP_pre_TSS.csv')

# Graph 2
colors <- colorRampPalette(c('white','red','gray','black'))(100)
layout(mat=matrix(c(1,2,0,3), 2, 2),
       widths=c(2,2,2),
       heights=c(0.5,5,0.5,5), TRUE)
par(mar=c(2.5,1,1.5,1))
image(seq(0, max(HPafterTSS), length.out=100), 1,
      matrix(seq(0, max(HPafterTSS), length.out=100),100,1),
      col = colors,
      xlab='Distance from TSS', ylab='',
      main='Number of reads', yaxt='n',
      lwd=3, axes=TRUE)
box(col='black', lwd=2)
image(x=seq(-1000, 1000, length.out=20),
      y=1:nrow(HPafterTSS),
      z=t(HPafterTSS[order(rowSums(HPafterTSS)),]),
      col=colors,
      xlab='Distance from TSS (bp)',
      ylab='Promoters', lwd=2)
box(col='black', lwd=2)
abline(v=0, lwd=1, col='gray')

plot(x=seq(-1000, 1000, length.out=20),
     y=colMeans(HPafterTSS),
     ty='b', pch=19,
     col='red4',lwd=2,
     ylab='Mean tag count',
     xlab='Distance from TSS (bp)')
abline(h=seq(1,100,by=5),
       v=seq(-1000, 1000, length.out=20),
       lwd=0.25, col='gray')
box(col='black', lwd=2)

# Graph 3
colors <- colorRampPalette(c('white','red','gray','black'))(100)
layout(mat=matrix(c(1,2,0,3), 2, 2),
       widths=c(2,2,2),
       heights=c(0.5,5,0.5,5), TRUE)
par(mar=c(2.5,1,1.5,1))
image(seq(0, max(HPpreTSS), length.out=100), 1,
      matrix(seq(0, max(HPpreTSS), length.out=100),100,1),
      col = colors,
      xlab='Distance from TSS', ylab='',
      main='Number of reads', yaxt='n',
      lwd=3, axes=TRUE)
box(col='black', lwd=2)
image(x=seq(-1000, 1000, length.out=20),
      y=1:nrow(HPpreTSS),
      z=t(HPpreTSS[order(rowSums(HPpreTSS)),]),
      col=colors,
      xlab='Distance from TSS (bp)',
      ylab='Promoters', lwd=2)
box(col='black', lwd=2)
abline(v=0, lwd=1, col='gray')

plot(x=seq(-1000, 1000, length.out=20),
     y=colMeans(HPpreTSS),
     ty='b', pch=19,
     col='red4',lwd=2,
     ylab='Mean tag count',
     xlab='Distance from TSS (bp)')
abline(h=seq(1,100,by=5),
       v=seq(-1000, 1000, length.out=20),
       lwd=0.25, col='gray')
box(col='black', lwd=2)

# seqrch genes by go
pre_genes <- go.table[go.table$ID %in% genespreTSS,]
selected_pre_genes <- pre_genes[grepl("cell cycle", pre_genes$GOTERM_BP_DIRECT)
                                | grepl("cell cycle", pre_genes$GOTERM_MF_DIRECT),]
selected_pre_genes <- egs[egs$ensembl_gene_id %in% selected_pre_genes$ID,] <- egs[egs$ensembl_gene_id %in% selected_pre_genes$ID,]
selected_pre_genes <- selected_pre_genes[!duplicated(selected_pre_genes$ensembl_gene_id),]
View(selected_pre_genes)