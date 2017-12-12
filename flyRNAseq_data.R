source("https://bioconductor.org/biocLite.R")
require(SummarizedExperiment)
require(DESeq2)
require(magrittr)
library(org.Dm.eg.db)
library(edgeR)
library(GeneOverlap)
library(psych)

setwd('/Volumes/NGSDATA/FlyRNAseq')
load('count.fly3.Rdata')

raw_counts <- count.fly3[-(1:2)] 
rownames(raw_counts) <- count.fly3[,1]

# Remove low expressed genes
raw_counts[is.na(raw_counts)] <- 0
myCPM <- cpm(raw_counts)
thresh <- myCPM > 20
keep <- rowSums(thresh) >= 2
# plot(myCPM[,1],raw_counts[,1],ylim=c(0,50),xlim=c(0,3))
# abline(h=10)
counts.keep <- raw_counts[keep,]

raw_counts_cad <- counts.keep["FBgn0003391",]
cor_mat <- cor(t(counts.keep),as.numeric(raw_counts_cad), use='pairwise.complete.obs')
cor_pvl <- corr.p(cor_mat , nrow(cor_mat), adjust="fdr", alpha=.05)

cor_mat <- cbind(cor_pvl$r, cor_pvl$p, cor_pvl$t)

a <- subset(cor_mat, cor_mat[,1] > 0.8)
b <- subset(cor_mat, cor_mat[,1] < -0.6)
all <- rbind(a,b)

all <- all[order(all[,1],decreasing=TRUE),]
all <- all[-1,]

names <- cbind(as.vector(count.fly3[,1]),as.vector(count.fly3[,2]))
rownames(names) <- count.fly3[,1]

all2 <- cbind(t(t(names[rownames(all),1:2])),all)
all2 <- all2[,-1]
write.csv(all2,'Cad_correlation.csv')

HP <- read.csv('genes0-200.csv')
Common <- t(t(intersect(HP[,6],noquote(rownames(all2)))))
rownames(Common) <- Common
Common <- all2[rownames(Common),]
write.csv(Common,'IntersectExpHP1.csv')

go.obj <- newGeneOverlap(HP[,6],noquote(rownames(all2)), genome.size=17918)
go.obj <- testGeneOverlap(go.obj)
