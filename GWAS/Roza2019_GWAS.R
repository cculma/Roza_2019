rm(list = ls())

library(GWASpoly)
library(tidyverse)
library(vcfR)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(data.table)
library(ggthemes)
library(hrbrthemes)
library(plotly)

####################
# GWASpoly using all mr.bean data (pheno_nph.csv) and genotypic matrix (MPP_Ms2_GWASPoly.txt)
setwd("~/Documents/git/big_files/")

pheno <- read.csv("BLUP_Yi_sqrt_SpATS_DArT.csv", row.names = 1)

pheno <- read.csv("BLUE_Yi_sqrt_SpATS_DArT.csv", row.names = 1)
head(pheno)
trait1 <- colnames(pheno)[1:(length(colnames(pheno))-3)]
trait1
params <- set.params(fixed=c("PC1","PC2","PC3"),
                     fixed.type=rep("numeric",3), n.PC = 3)
models_1 <- c("general", "additive", "1-dom", "2-dom",  "diplo-additive", "diplo-general")

# "ACGT"
# numeric
# Roza2019_06_GWASPoly.txt # DArTag

data_1 <- read.GWASpoly(ploidy=4, 
                        pheno.file="BLUE_Yi_sqrt_SpATS_DArT.csv", 
                        geno.file="Roza2019_06_GWASPoly.txt", 
                        format="numeric", n.traits=length(trait1), delim=",")

data_2 <- set.K(data = data_1, LOCO = T, n.core = 16)
Yi_data_3 <- GWASpoly(data = data_2, models = models_1, traits = trait1, params = params, n.core = 16)

save(Yi_data_3, file = "~/Documents/git/big_files/Yi_st2_51081_log.RData")
save(Yi_data_3, file = "~/Documents/git/big_files/Yi_st1_51081_sqrt.RData")


save(Yi_data_3, file = "~/Documents/git/big_files/Yi_st1_51081.RData")

# load("~/Documents/git/big_files/Yi_st1_51081.RData")
# load("~/Documents/git/big_files/Yi_st2_51081.RData")

data_5 <- set.threshold(Yi_data_3, method= "Bonferroni", level=0.05)
QTL_01 <- get.QTL(data_5)
# manhattan.plot(data = data_5)

QTL_02 <- QTL_01 %>% distinct(QTL_01$Marker, .keep_all = T)

lev0 <- unique(QTL_01$Trait)

lev0 <- subset(QTL_02$Trait, grepl("ST0_", QTL_02$Trait))
lev1 <- subset(QTL_02$Trait, grepl("ST1_", QTL_02$Trait))

lev4 <- subset(QTL_01$Trait, grepl("ST4_", QTL_01$Trait))

# save.image("~/Documents/Cesar/blup_data/Roza2019/Analysis_2021/GWAS/data_3_80177_year.RData")
# load("~/Documents/Cesar/blup_data/Roza2019/Analysis_2021/GWAS/data_3_80177_year.RData")

M0 <- manhattan.plot(data = data_5, traits = lev0) + theme_classic(base_family = "Arial", base_size = 12) + theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y = element_text(size = 12), plot.tag = element_blank()) 

ggsave(filename = "~/Documents/git/big_files/M0_sqrt.jpg", plot = M0, width = 16, height = 16)


M1 <- manhattan.plot(data = data_5, traits = lev0) + theme_classic(base_family = "Arial", base_size = 12) + theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y = element_text(size = 12), plot.tag = element_blank()) 

ggsave(filename = "~/Documents/git/big_files/M1_sqrt.jpg", plot = M1, width = 16, height = 16)


manhattan.plot(data = data_5, traits = "ST4_Yi") + theme_classic(base_family = "Arial", base_size = 12) + theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y = element_text(size = 12), plot.tag = element_blank())
