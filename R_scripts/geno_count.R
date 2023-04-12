# This script allow to count markers by chromosome

library(tidyverse)
geno <- read.csv("~/OneDrive - Washington State University (email.wsu.edu)/Roza_2019/GWASpoly_results/MPP_Ms2_GWASPoly.txt", row.names = 1)
geno[1:5, 1:5]
geno[82150:82156, 1:5]
dim(geno)
geno$Chrom[grepl("^contig*", geno$Chrom)] <- "Contig"
cc <- count(geno, Chrom)
sum(cc$n)

# This script allow to summarize redundant associated markers by GWASpoly model and by Chromosome

load("~/OneDrive - Washington State University (email.wsu.edu)/Roza_2019/GWASpoly_results/data_3_80177_year.RData")

QTL_02 <- QTL_01 %>% filter(Trait != "MET_all")
QTL_02$Chrom <- as.character(QTL_02$Chrom)
QTL_02$Chrom[grepl("^contig*", QTL_02$Chrom)] <- "Contig"
df1 <- count(QTL_02, Model)
df2 <- count(QTL_02, Chrom)

