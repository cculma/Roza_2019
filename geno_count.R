# count markers by chromosome
library(tidyverse)
geno <- read.csv("~/OneDrive - Washington State University (email.wsu.edu)/Roza_2019/GWASpoly_results/MPP_Ms2_GWASPoly.txt", row.names = 1)
geno[1:5, 1:5]
geno[82150:82156, 1:5]
dim(geno)
geno$Chrom[grepl("^contig*", geno$Chrom)] <- "Contig"
cc <- count(geno, Chrom)
sum(cc$n)