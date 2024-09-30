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

pheno <- read.csv("BLUP_Yi_sqrt_SpATS_DArT.csv", row.names = 1) # ST2

head(pheno)
trait1 <- colnames(pheno)[1:(length(colnames(pheno))-3)]
trait1
# overall by month: # "yi_may"  "yi_jun"  "yi_jul"  "yi_aug"  "yi_sep"  
# overall by year: "yi_2020" "yi_2021" "yi_2022" "yi_2023"
# overall all months and all years: "yi_hrv"  "yi_year"

params <- set.params(fixed=c("PC1","PC2","PC3"),
                     fixed.type=rep("numeric",3), n.PC = 3)
models_1 <- c("general", "additive", "1-dom", "2-dom",  "diplo-additive", "diplo-general")

# numeric
# Roza2019_06_GWASPoly.txt # DArTag

data_1 <- read.GWASpoly(ploidy=4, 
                        pheno.file="BLUP_Yi_sqrt_SpATS_DArT.csv", 
                        geno.file="Roza2019_06_GWASPoly.txt", 
                        format="numeric", n.traits=length(trait1), delim=",")

data_2 <- set.K(data = data_1, LOCO = T, n.core = 16)
Yi_data_3 <- GWASpoly(data = data_2, models = models_1, traits = trait1, params = params, n.core = 16)

save(Yi_data_3, file = "~/Documents/git/big_files/Yi_st2_51081_sqrt.RData") # data sqrt transformed yield

# load("~/Documents/git/big_files/Yi_st2_51081_sqrt.RData")

data_5 <- set.threshold(Yi_data_3, method= "Bonferroni", level=0.05)
QTL_01 <- get.QTL(data_5)

QTL_02 <- QTL_01 %>% distinct(QTL_01$Marker, .keep_all = T)

lev0 <- unique(QTL_01$Trait)

M1 <- manhattan.plot(data = data_5, traits = lev0) + theme_classic(base_family = "Arial", base_size = 12) + theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y = element_text(size = 12), plot.tag = element_blank()) 
ggsave(filename = "~/Documents/git/big_files/M1_sqrt.jpg", plot = M1, width = 16, height = 16)


cc <- dplyr::count(QTL_01, Trait)

QTL_02 <- QTL_01 %>% dplyr::filter(Trait == "yi_may")
QTL_02 <- QTL_01 %>% dplyr::filter(Trait == "yi_jun")
QTL_02 <- QTL_01 %>% dplyr::filter(Trait == "yi_jul")
QTL_02 <- QTL_01 %>% dplyr::filter(Trait == "yi_aug")
QTL_02 <- QTL_01 %>% dplyr::filter(Trait == "yi_sep")

QTL_02 <- QTL_01 %>% dplyr::filter(Trait == "yi_2020")
QTL_02 <- QTL_01 %>% dplyr::filter(Trait == "yi_2021")
QTL_02 <- QTL_01 %>% dplyr::filter(Trait == "yi_2022")
QTL_02 <- QTL_01 %>% dplyr::filter(Trait == "yi_2023")

QTL_02 <- QTL_01 %>% dplyr::filter(Trait == "yi_hrv")
QTL_02 <- QTL_02 %>% distinct(Marker, .keep_all = T)
