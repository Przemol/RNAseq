
printSig <- function(R, x=0.6) {
  library(tibble)
  library(dplyr)
  library(DESeq2)
  
  if(x>0){
    sig <- R[R>x,]
    tbl <- rev(sort(sig)) %>% as.data.frame() %>%  rownames_to_column('entrez') %>% tbl_df 
  } else {
    sig <- R[R<x,]
    tbl <- sort(sig) %>% as.data.frame() %>%  rownames_to_column('entrez') %>% tbl_df 
  }
  colnames(tbl) = c('entrez', 'cor')
  
  tbl$SYMBOL <- suppressMessages( mapIds(org.Hs.eg.db, tbl$entrez, "SYMBOL", "ENTREZID") )
  tbl$GENENAME <- suppressMessages( mapIds(org.Hs.eg.db, tbl$entrez, "GENENAME", "ENTREZID") )
  return(tbl)
}



getHiCor <- function(SE, SYMBOL='CDH1', rpkm='TRUE', cor_cutoff=0.6, anticor_cutoff=-0.6) {
  
  ENTREZID <- suppressMessages( mapIds(org.Hs.eg.db, SYMBOL, "ENTREZID", "SYMBOL") )

  message('Assesing ', SYMBOL, ' (ENTREZID: ', ENTREZID, ')')
  
  if(rpkm) {
    M <- fpkm(DESeqDataSet(SE, ~1), robust = TRUE)
  } else {
    M <- assay(SE)
  }

  Mcand <- M[ENTREZID,]
  cor_mat <- cor(t(M), Mcand, use='pairwise.complete.obs')
  
  a <- printSig(cor_mat, cor_cutoff)
  b <- printSig(cor_mat, anticor_cutoff)
  
  sig <- rbind(a,b)
  
  return(sig)
}


