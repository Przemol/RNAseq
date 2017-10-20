plotNet <- function(data=NULL, cc=NULL, method="spearman", lasso=TRUE, rmD=TRUE, prefix='Genes and miRNA:', ...) {
  library(qgraph)
  
  if(is.null(cc)) {
    cc <- cor(data, use = 'pairwise.complete.obs', method=method)
  }
  title <- paste(prefix, method, 'correlation coefiicient')
  
  if(lasso) {
    pd <- matrixcalc::is.positive.definite(cc)
    if(!pd) {
      warning('Matrix is not PD, finding nearPD')
      cc <- nearPD(cc, corr = TRUE)$mat %>% as.matrix
    }
    EBICgraph <- EBICglasso(cc, 20000, 0.5)
    title <- paste(title, 'with EBICglasso')
  } else {
    if(rmD) diag(cc) <- NA
    EBICgraph <- cc
  }
  
  z <- qgraph(EBICgraph, layout = "spring", title = title, DoNotPlot=TRUE, ...)
  z$graphAttributes$Nodes$labels <- names(z$graphAttributes$Nodes$labels)
  plot(z)
}

plotCC <- function(data, cutoff=1, method="spearman") {
  library(corrplot)
  
  cc <- cor(data, use = 'pairwise.complete.obs', method=method)
  corrplot(cc, order = "AOE", type = 'lower', diag=T, method = 'number', tl.pos = 'lt', p.mat=(cc >= cutoff & !diag(ncol(cc))))
  corrplot(cc, order = "AOE", type = 'upper',diag=F, add = TRUE, method = 'circle', tl.pos = 'n', p.mat=(cc >= cutoff))
}

plotCor <- function(CorMat) {
  library(reshape2)
  library(GGally)
  
  reorder_cormat <- function(cormat){
    # Use correlation between variables as distance
    dd <- as.dist((1-cormat)/2)
    hc <- hclust(dd)
    cormat <-cormat[hc$order, hc$order]
  }
  
  CC <- reorder_cormat(CorMat)
  #CC <- reorder_cormat(MM)
  #CC <- reorder_cormat(EBICgraph)
  
  
  hmap <- ggcorr(
    data=NULL, cor_matrix = as.data.frame(CC[rev(rownames(CC)),rev(rownames(CC))]),
    limits = FALSE, hjust = 0.7, vjust = .7, label = TRUE, layout.exp=1,
    label_round = 2, angle = 0, label_size = 5, legend.size = 12, size=5.5,
    low = "darkred", mid = "white", high = "darkgreen"
    #,mapping=aes(c("H3K9me2", "H3K9me3", "HPL-2", "MET-2", "LET-418", "LIN-13", "LIN-61"))
  )
  
  
  hmap <- hmap + theme(
    legend.justification = c(1, 0),
    legend.position = c(0.4, 0.6),
    legend.direction = "horizontal",
    legend.title = element_text(size = 14)
  ) + guides(fill = guide_colorbar(
    title="Correlation",
    barwidth = 7, barheight = 1,
    title.position = "top",
    title.hjust = 0.5)
  );
  return(hmap)
}

collapsemaxcor <- function(MAT, method="spearman") {
  CorMat_p <- cor(MAT, use = 'pairwise.complete.obs', method=method)
  diag(CorMat_p) <- 0
  ind <- rownames(which(CorMat_p==max(CorMat_p), arr.ind = TRUE))
  
  message('Collapsing: ', paste(ind, collapse = " and "), 'at cc=', max(CorMat_p))
  
  MM <- MAT[, !colnames(MAT) %in% ind]
  rs <- rowSums(MAT[, ind], na.rm = TRUE)
  rs[rs==0] <- NA
  MM <- cbind(MM, rs)
  colnames(MM)[ncol(MM)] <- paste(ind, collapse = '\n')
  
  PD <- matrixcalc::is.positive.definite(cor(MM, use = 'pairwise.complete.obs', method=method))
  if(!PD) message('Still not PD')
  attr(MM, 'PD') <- PD
  attr(MM, 'crm') <- max(CorMat_p)
  return(MM)
}


