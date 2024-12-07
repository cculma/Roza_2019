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
setwd("~/Documents/git/big_files/")

pheno <- read.csv("BLUE_Yi_sqrt_SpATS_DArT.csv", row.names = 1) # ST1 (SpATS)

head(pheno)
trait1 <- colnames(pheno)[1:(length(colnames(pheno))-3)]
trait1
params <- set.params(fixed=c("PC1","PC2","PC3"),
                     fixed.type=rep("numeric",3), n.PC = 3)
models_1 <- c("general", "additive", "1-dom", "2-dom",  "diplo-additive", "diplo-general")

# numeric
# Roza2019_06_GWASPoly.txt # DArTag

data_1 <- read.GWASpoly(ploidy=4, 
                        pheno.file="BLUE_Yi_sqrt_SpATS_DArT.csv", 
                        geno.file="Roza2019_06_GWASPoly.txt", 
                        format="numeric", n.traits=length(trait1), delim=",")

data_2 <- set.K(data = data_1, LOCO = T, n.core = 16)
Yi_data_3 <- GWASpoly(data = data_2, models = models_1, traits = trait1, params = params, n.core = 16)

# save(Yi_data_3, file = "~/Documents/git/big_files/Yi_st1_51081_sqrt.RData") # data sqrt transformed

# load("~/Documents/git/big_files/Yi_st1_51081_sqrt.RData")

data_5 <- set.threshold(Yi_data_3, method= "Bonferroni", level=0.05)
QTL_01 <- get.QTL(data_5)



# lev0 <- unique(QTL_01$Trait)

M0 <- manhattan.plot(data = data_5, traits = lev0) + theme_classic(base_family = "Arial", base_size = 12) + theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y = element_text(size = 12), plot.tag = element_blank()) 
ggsave(filename = "~/Documents/git/big_files/M0_sqrt.jpg", plot = M0, width = 16, height = 16)



# count markers -----------------------------------------------------------

cc <- dplyr::count(QTL_01, Trait)


QTL_02 <- QTL_01 %>% dplyr::filter(Trait == "may_20")
QTL_02 <- QTL_01 %>% dplyr::filter(Trait == "may_21")
QTL_02 <- QTL_01 %>% dplyr::filter(Trait == "may_22")
QTL_02 <- QTL_01 %>% dplyr::filter(Trait == "may_23")

QTL_02 <- QTL_01 %>% dplyr::filter(Trait == "jun_20")
QTL_02 <- QTL_01 %>% dplyr::filter(Trait == "jun_21")
QTL_02 <- QTL_01 %>% dplyr::filter(Trait == "jun_22")
QTL_02 <- QTL_01 %>% dplyr::filter(Trait == "jun_23")

QTL_02 <- QTL_01 %>% dplyr::filter(Trait == "jul_20")
QTL_02 <- QTL_01 %>% dplyr::filter(Trait == "jul_21")
QTL_02 <- QTL_01 %>% dplyr::filter(Trait == "jul_22")
QTL_02 <- QTL_01 %>% dplyr::filter(Trait == "jul_23")

QTL_02 <- QTL_01 %>% dplyr::filter(Trait == "aug_20")
QTL_02 <- QTL_01 %>% dplyr::filter(Trait == "aug_21")
QTL_02 <- QTL_01 %>% dplyr::filter(Trait == "aug_22")
QTL_02 <- QTL_01 %>% dplyr::filter(Trait == "aug_23")

QTL_02 <- QTL_01 %>% dplyr::filter(Trait == "sep_20")
QTL_02 <- QTL_01 %>% dplyr::filter(Trait == "sep_21")
QTL_02 <- QTL_01 %>% dplyr::filter(Trait == "sep_22")
QTL_02 <- QTL_01 %>% dplyr::filter(Trait == "sep_23")
QTL_02 <- QTL_02 %>% distinct(Marker, .keep_all = T)


