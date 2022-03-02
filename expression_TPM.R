# expression levels
# load TPM matrix using file cpms.3.2 generated with the script tximport.R in iso_seq_shen repository
rm(list = ls())

library(tidyverse)
library(pheatmap)
library(DESeq2)


load("~/OneDrive - Washington State University (email.wsu.edu)/Roza_2019/GWASpoly_results/data_3_80177_year.RData")
load("~/OneDrive - Washington State University (email.wsu.edu)/RNA/Ms_Shen/cpms.3.2.RData")
load("~/OneDrive - Washington State University (email.wsu.edu)/RNA/Ms_Shen/QTL_10.RData")
load("~/OneDrive - Washington State University (email.wsu.edu)/RNA/Ms_Shen/tximport.RData")

lev_2 <- c("G25710", ) 

QTL_09 <- QTL_08 %>% group_by(gene_id) %>% summarise(Marker1 = paste(Marker1, collapse = ";")) 
QTL_10 <- QTL_08 %>% distinct(gene_id, .keep_all = TRUE) %>% dplyr::select(6:8)
QTL_10 <- inner_join(QTL_09, QTL_10, by = "gene_id")
save(QTL_10, file = "~/OneDrive - Washington State University (email.wsu.edu)/RNA/Ms_Shen/QTL_10.RData")

QTL_10$gene_id
colnames(cpms.3.1)
colnames(cpms.3.2)

# Sanarac DS
cpms.3.3 <- cpms.3.2[,c(1:9)]
colnames(cpms.3.3)

cpms.4 <- left_join(QTL_10, cpms.3.2, by = "gene_id") %>% dplyr::select(5,7:18) %>% remove_rownames() %>% column_to_rownames(var = "lid") %>% filter_all(any_vars(. != 0))
colnames(cpms.4)

cpms.5 <- cpms.4 %>% tibble::rownames_to_column(var = "lid") %>% separate(1, c("gene_id", "isoform"), sep = "\\.", remove = FALSE, convert = FALSE, extra = "warn") %>% group_by(gene_id) %>% summarise_at(vars(4:13), mean) %>% remove_rownames() %>% column_to_rownames(var = "gene_id")

colnames(cpms.5)


# cc <- dplyr::count(cpms.5, gene_id)
# cc1 <- cc %>% dplyr::filter(n >= 100)

# cpms.6 <- left_join(cc1, cpms.3.3, by = "gene_id") %>% dplyr::select(3,5:10) %>% remove_rownames() %>% column_to_rownames(var = "lid") %>% filter_all(any_vars(. != 0))
colnames(cpms.6)


pheatmap(cpms.5, drop_levels = T, show_rownames = F)

cal_z_score <- function(x){(x - mean(x)) / sd(x)}
data_subset_norm <- t(apply(cpms.5, 1, cal_z_score))
pheatmap(data_subset_norm, drop_levels = T, show_rownames = F)

data_subset_norm[1:5,1:5]
dim(data_subset_norm)


ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~ condition)


# to remove rows where all values are 0
# cpms.4 <- cpms.4[rowSums(cpms.4[]) > 0,]

# to normalize data in log value
# cpms.7 <- log2(cpms.4 + 1)
colnames(cpms.4)
cpms.7 <- cpms.4[,c(10,7,4,1)]
data_subset_norm <- t(apply(cpms.7, 1, cal_z_score))
data_subset_norm <- na.omit(data_subset_norm)
pheatmap(data_subset_norm, drop_levels = T, show_rownames = F)
