rm(list = ls()) # clean Global Environment
# setwd("~/Documents/Cesar/blup_data/Roza2019/Analysis_2021/GWAS/")
# setwd("~/OneDrive - Washington State University (email.wsu.edu)/Roza_2019/git/Roza2019/")

library(GWASpoly)
library(tidyverse)
library(vcfR)
library(parallel)
library(doParallel)
library(iterators)
library(foreach)
library(tidyr)
library(devtools)
library(sommer)
library(ggplot2)
library(ggpubr)
library(data.table)
library(ggthemes)
library(hrbrthemes)
library(VennDiagram)
library(plotly)

# load("~/OneDrive - Washington State University (email.wsu.edu)/Roza_2019/GWASpoly_results/data_3_80177_year.RData")
# load("~/Documents/Cesar/blup_data/Roza2019/Analysis_2021/GWAS/data_3_80177_year.RData")

####################
# GWASpoly using all mr.bean data (pheno_nph.csv) and genotypic matrix (MPP_Ms2_GWASPoly.txt)

setwd("~/Documents/Cesar/blup_data/Roza2019/Analysis_2021/GWAS/")
pheno <- read.csv("~/Documents/Cesar/blup_data/Roza2019/Analysis_2021/GWAS/pheno_yi.csv", row.names = 1)
head(pheno)
trait1 <- colnames(pheno)[1:(length(colnames(pheno))-3)]
trait1
params <- set.params(fixed=c("PC1","PC2","PC3"),
                     fixed.type=rep("numeric",3), n.PC = 3)
models_1 <- c("general", "additive", "1-dom", "2-dom",  "diplo-additive", "diplo-general")
# models_1 <- c("general", "additive", "1-dom",  "diplo-additive")
data_1 <- read.GWASpoly(ploidy=4, 
                        pheno.file="pheno_yi.csv", 
                        geno.file="MPP_Ms2_GWASPoly.txt", 
                        format="ACGT", n.traits=length(trait1), delim=",")

data_2 <- set.K(data = data_1, LOCO = T, n.core = 32)
Yi_data_3 <- GWASpoly(data = data_2, models = models_1, traits = trait1, params = params, n.core = 32)
save(Yi_data_3, file = "~/Documents/Cesar/blup_data/Roza2019/Analysis_2021/GWAS/Yi_data_3_82156.RData")

data_5 <- set.threshold(data_3, method= "Bonferroni", level=0.05)
QTL_01 <- get.QTL(data_5)
QTL_02 <- QTL_01 %>% distinct(Marker, .keep_all = T) 

# save.image("~/Documents/Cesar/blup_data/Roza2019/Analysis_2021/GWAS/data_3_80177_year.RData")
# load("~/Documents/Cesar/blup_data/Roza2019/Analysis_2021/GWAS/data_3_80177_year.RData")
