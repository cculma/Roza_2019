# expression levels
# load TPM matrix using file cpms.3.2 generated with the script tximport.R in iso_seq_shen repository
rm(list = ls())

library(tidyverse)
library(pheatmap)

load("~/OneDrive - Washington State University (email.wsu.edu)/Roza_2019/GWASpoly_results/data_3_80177_year.RData")
load("~/OneDrive - Washington State University (email.wsu.edu)/RNA/Ms_Shen/cpms.3.2.RData")
load("~/OneDrive - Washington State University (email.wsu.edu)/RNA/Ms_Shen/QTL_10.RData")
load("~/OneDrive - Washington State University (email.wsu.edu)/RNA/Ms_Shen/tximport.RData")



QTL_09 <- QTL_08 %>% group_by(gene_id) %>% summarise(Marker1 = paste(Marker1, collapse = ";")) 
QTL_10 <- QTL_08 %>% distinct(gene_id, .keep_all = TRUE) %>% dplyr::select(6:8)
QTL_10 <- inner_join(QTL_09, QTL_10, by = "gene_id")
save(QTL_10, file = "~/OneDrive - Washington State University (email.wsu.edu)/RNA/Ms_Shen/QTL_10.RData")

colnames(cpms.4)
cpms.4 <- left_join(QTL_10, cpms.3.2, by = "gene_id") %>% dplyr::select(5,7:18) %>% remove_rownames() %>% column_to_rownames(var = "lid")

cpms.4 <- cpms.4[rowSums(cpms.4[])>0,]
cpms.5 <- log2(cpms.4 + 1)
cpms.6 <- cpms.5[apply(cpms.5, 1, var) > 0, ]
dim(cpms.4)
dim(cpms.6)

pheatmap(cpms.4)

cal_z_score <- function(x){(x - mean(x)) / sd(x)}
data_subset_norm <- t(apply(cpms.6, 1, cal_z_score))
pheatmap(data_subset_norm, drop_levels = T, show_rownames = F)

data_subset_norm[1:5,1:5]
dim(data_subset_norm)
