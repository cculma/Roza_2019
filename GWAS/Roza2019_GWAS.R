rm(list = ls())
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
setwd("~/Documents/blup_data/Roza2019/Analysis_2021/GWAS/")

pheno <- read.csv("pheno_yi.csv", row.names = 1)
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
Yi_data_3 <- GWASpoly(data = data_2, models = models_1, traits = trait1, params = params, n.core = 32)
# save(Yi_data_3, file = "~/Documents/blup_data/Roza2019/Analysis_2021/GWAS/Yi_data_3_82156.RData")
load("~/Documents/blup_data/Roza2019/Analysis_2021/GWAS/Yi_data_3_82156.RData")

data_5 <- set.threshold(Yi_data_3, method= "Bonferroni", level=0.05)
QTL_01 <- get.QTL(data_5)
QTL_02 <- QTL_01 %>% distinct(Trait, .keep_all = T) 
QTL_02$Trait

QTL_04 <- QTL_01 %>% dplyr::filter(Trait == "ST4_Yi")

lev0 <- subset(QTL_02$Trait, grepl("ST0_", QTL_02$Trait))
lev1 <- subset(QTL_02$Trait, grepl("ST1_", QTL_02$Trait))

lev4 <- subset(QTL_01$Trait, grepl("ST4_", QTL_01$Trait))

# save.image("~/Documents/Cesar/blup_data/Roza2019/Analysis_2021/GWAS/data_3_80177_year.RData")
# load("~/Documents/Cesar/blup_data/Roza2019/Analysis_2021/GWAS/data_3_80177_year.RData")

M0 <- manhattan.plot(data = data_5, traits = lev0) + theme_classic(base_family = "Arial", base_size = 12) + theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y = element_text(size = 12), plot.tag = element_blank()) 

ggsave(filename = "~/Documents/Cesar/blup_data/Roza2019/Analysis_2021/GWAS/M0.jpg", plot = M0, width = 16, height = 16)

M1 <- manhattan.plot(data = data_5, traits = lev1) + theme_classic(base_family = "Arial", base_size = 12) + theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y = element_text(size = 12), plot.tag = element_blank()) 

ggsave(filename = "~/Documents/Cesar/blup_data/Roza2019/Analysis_2021/GWAS/M1.jpg", plot = M1, width = 16, height = 16)



manhattan.plot(data = data_5, traits = "ST4_Yi") + theme_classic(base_family = "Arial", base_size = 12) + theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y = element_text(size = 12), plot.tag = element_blank())
