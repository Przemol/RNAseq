counts.fly <- read.table('gene_rpkm_report_fb_2017_05.tsv', sep = '\t', header = FALSE, fill = TRUE)

install.packages("dplyr")
library("dplyr")
count.fly2 <- select(counts.fly, V2, V3, V7, V8)

count.fly2 <- split(count.fly2, count.fly2$V7, drop=F)

count.fly3 <- merge(count.fly2$`mE_mRNA_em0-2hr`, count.fly2$`mE_mRNA_em2-4hr`, by="V2", all = T)
count.fly3 <- select(count.fly3, V2, V3.x, V8.x, V8.y) 
colnames(count.fly3) <- c('V2', 'Name','mE_mRNA_em0-2hr', 'mE_mRNA_em2-4hr') -> cn

count.fly2 <- select(counts.fly, V2, V7, V8)
count.fly2 <- split(count.fly2, count.fly2$V7, drop=F)

count.fly3 <- merge(count.fly3, count.fly2$`mE_mRNA_em4-6hr`, by="V2", all = T)
count.fly3 <- count.fly3[,-5]
colnames(count.fly3) <- c(cn, 'mE_mRNA_em4-6hr') -> cn

count.fly3 <- merge(count.fly3, count.fly2$`mE_mRNA_em6-8hr`, by="V2", all = T)
count.fly3 <- count.fly3[,-6]
colnames(count.fly3) <- c(cn, 'mE_mRNA_em6-8hr') -> cn

count.fly3 <- merge(count.fly3, count.fly2$`mE_mRNA_em8-10hr`, by="V2", all = T)
count.fly3 <- count.fly3[,-7]
colnames(count.fly3) <- c(cn, 'mE_mRNA_em8-10hr') -> cn

count.fly3 <- merge(count.fly3, count.fly2$`mE_mRNA_em10-12hr`, by="V2", all = T)
count.fly3 <- count.fly3[,-8]
colnames(count.fly3) <- c(cn, 'mE_mRNA_em10-12hr') -> cn

count.fly3 <- merge(count.fly3, count.fly2$`mE_mRNA_em12-14hr`, by="V2", all = T)
count.fly3 <- count.fly3[,-9]
colnames(count.fly3) <- c(cn, 'mE_mRNA_em12-14hr') -> cn

count.fly3 <- merge(count.fly3, count.fly2$`mE_mRNA_em14-16hr`, by="V2", all = T)
count.fly3 <- count.fly3[,-10]
colnames(count.fly3) <- c(cn, 'mE_mRNA_em14-16hr') -> cn

count.fly3 <- merge(count.fly3, count.fly2$`mE_mRNA_em16-18hr`, by="V2", all = T)
count.fly3 <- count.fly3[,-11]
colnames(count.fly3) <- c(cn, 'mE_mRNA_em16-18hr') -> cn

count.fly3 <- merge(count.fly3, count.fly2$`mE_mRNA_em18-20hr`, by="V2", all = T)
count.fly3 <- count.fly3[,-12]
colnames(count.fly3) <- c(cn, 'mE_mRNA_em18-20hr') -> cn

count.fly3 <- merge(count.fly3, count.fly2$`mE_mRNA_em20-22hr`, by="V2", all = T)
count.fly3 <- count.fly3[,-13]
colnames(count.fly3) <- c(cn, 'mE_mRNA_em20-22hr') -> cn

count.fly3 <- merge(count.fly3, count.fly2$`mE_mRNA_em22-24hr`, by="V2", all = T)
count.fly3 <- count.fly3[,-14]
colnames(count.fly3) <- c(cn, 'mE_mRNA_em22-24hr') -> cn

count.fly3 <- merge(count.fly3, count.fly2$`mE_mRNA_L1`, by="V2", all = T)
count.fly3 <- count.fly3[,-15]
colnames(count.fly3) <- c(cn, 'mE_mRNA_L1') -> cn

count.fly3 <- merge(count.fly3, count.fly2$`mE_mRNA_L2`, by="V2", all = T)
count.fly3 <- count.fly3[,-16]
colnames(count.fly3) <- c(cn, 'mE_mRNA_L2') -> cn

count.fly3 <- merge(count.fly3, count.fly2$`mE_mRNA_L3_12hr`, by="V2", all = T)
count.fly3 <- count.fly3[,-17]
colnames(count.fly3) <- c(cn, 'mE_mRNA_L3_12hr') -> cn

count.fly3 <- merge(count.fly3, count.fly2$`mE_mRNA_L3_PS1-2`, by="V2", all = T)
count.fly3 <- count.fly3[,-18]
colnames(count.fly3) <- c(cn, 'mE_mRNA_L3_PS1-2') -> cn

count.fly3 <- merge(count.fly3, count.fly2$`mE_mRNA_L3_PS3-6`, by="V2", all = T)
count.fly3 <- count.fly3[,-19]
colnames(count.fly3) <- c(cn, 'mE_mRNA_L3_PS3-6') -> cn

count.fly3 <- merge(count.fly3, count.fly2$`mE_mRNA_L3_PS7-9`, by="V2", all = T)
count.fly3 <- count.fly3[,-20]
colnames(count.fly3) <- c(cn, 'mE_mRNA_L3_PS7-9') -> cn

count.fly3 <- merge(count.fly3, count.fly2$`mE_mRNA_WPP`, by="V2", all = T)
count.fly3 <- count.fly3[,-21]
colnames(count.fly3) <- c(cn, 'mE_mRNA_WPP') -> cn

count.fly3 <- merge(count.fly3, count.fly2$`mE_mRNA_A_MateF_4d_ovary`, by="V2", all = T)
count.fly3 <- count.fly3[,-22]
colnames(count.fly3) <- c(cn, 'mE_mRNA_A_MateF_4d_ovary') -> cn

count.fly3 <- merge(count.fly3, count.fly2$`mE_mRNA_A_MateF_4d_head`, by="V2", all = T)
count.fly3 <- count.fly3[,-23]
colnames(count.fly3) <- c(cn, 'mE_mRNA_A_MateF_4d_head') -> cn

count.fly3 <- merge(count.fly3, count.fly2$`mE_mRNA_L3_CNS`, by="V2", all = T)
count.fly3 <- count.fly3[,-24]
colnames(count.fly3) <- c(cn, 'mE_mRNA_L3_CNS') -> cn

count.fly3 <- merge(count.fly3, count.fly2$`mE_mRNA_L3_Wand_carcass`, by="V2", all = T)
count.fly3 <- count.fly3[,-25]
colnames(count.fly3) <- c(cn, 'mE_mRNA_L3_Wand_carcass') -> cn

count.fly3 <- merge(count.fly3, count.fly2$`mE_mRNA_L3_Wand_dig_sys`, by="V2", all = T)
count.fly3 <- count.fly3[,-26]
colnames(count.fly3) <- c(cn, 'mE_mRNA_L3_Wand_dig_sys') -> cn

count.fly3 <- merge(count.fly3, count.fly2$`mE_mRNA_L3_Wand_fat`, by="V2", all = T)
count.fly3 <- count.fly3[,-27]
colnames(count.fly3) <- c(cn, 'mE_mRNA_L3_Wand_fat') -> cn

count.fly3 <- merge(count.fly3, count.fly2$`mE_mRNA_L3_Wand_imag_disc`, by="V2", all = T)
count.fly3 <- count.fly3[,-28]
colnames(count.fly3) <- c(cn, 'mE_mRNA_L3_Wand_imag_disc') -> cn

count.fly3 <- merge(count.fly3, count.fly2$`mE_mRNA_L3_Wand_saliv`, by="V2", all = T)
count.fly3 <- count.fly3[,-29]
colnames(count.fly3) <- c(cn, 'mE_mRNA_L3_Wand_saliv') -> cn

count.fly3 <- merge(count.fly3, count.fly2$`mE_mRNA_AdF_Ecl_1days`, by="V2", all = T)
count.fly3 <- count.fly3[,-30]
colnames(count.fly3) <- c(cn, 'mE_mRNA_AdF_Ecl_1days') -> cn

count.fly3 <- merge(count.fly3, count.fly2$`mE_mRNA_AdF_Ecl_5days`, by="V2", all = T)
count.fly3 <- count.fly3[,-31]
colnames(count.fly3) <- c(cn, 'mE_mRNA_AdF_Ecl_5days') -> cn

count.fly3 <- merge(count.fly3, count.fly2$`mE_mRNA_AdF_Ecl_30days`, by="V2", all = T)
count.fly3 <- count.fly3[,-32]
colnames(count.fly3) <- c(cn, 'mE_mRNA_AdF_Ecl_30days') -> cn

count.fly3 <- merge(count.fly3, count.fly2$`mE_mRNA_AdM_Ecl_1days`, by="V2", all = T)
count.fly3 <- count.fly3[,-33]
colnames(count.fly3) <- c(cn, 'mE_mRNA_AdM_Ecl_1days') -> cn

count.fly3 <- merge(count.fly3, count.fly2$`mE_mRNA_AdM_Ecl_5days`, by="V2", all = T)
count.fly3 <- count.fly3[,-34]
colnames(count.fly3) <- c(cn, 'mE_mRNA_AdM_Ecl_5days') -> cn

count.fly3 <- merge(count.fly3, count.fly2$`mE_mRNA_AdM_Ecl_30days`, by="V2", all = T)
count.fly3 <- count.fly3[,-35]
colnames(count.fly3) <- c(cn, 'mE_mRNA_AdM_Ecl_30days') -> cn

save(count.fly3, file='count.fly3.Rdata')