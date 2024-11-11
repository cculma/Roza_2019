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
library(GenomicRanges)
library(genomation)
library(plyranges)
library(Repitools)
library(devtools)
library(sommer)

####################
setwd("~/Documents/git/big_files/")

pheno <- read.csv("yield_DS.csv")
head(pheno)
trait1 <- c("predicted.value")
trait1
?set.params
params <- set.params(fixed=c("PC1","PC2","PC3"),
                     fixed.type= c("numeric","numeric","numeric"), n.PC = 3)

models_1 <- c("general", "additive", "1-dom", "2-dom",  "diplo-additive", "diplo-general")

data_1 <- read.GWASpoly(ploidy=4, 
                        pheno.file="yield_DS.csv", 
                        geno.file="Roza2019_06_GWASPoly.txt", 
                        format="numeric", n.traits=length(trait1), delim=",")


data_2 <- set.K(data = data_1, LOCO = T, n.core = 10)
Yi_data_3 <- GWASpoly(data = data_2, models = models_1, traits = trait1, params = params, n.core = 10)

save(Yi_data_3, file = "~/Documents/git/big_files/yield_DS.RData")
# save.image("~/Documents/git/big_files/yield_DS.RData")
# load("~/Documents/git/big_files/data_Yi_DS_20.RData")

data_5 <- set.threshold(Yi_data_3, method= "Bonferroni", level=0.05)
QTL_01 <- get.QTL(data_5)

QTL_02 <- QTL_01 %>% distinct(Marker, .keep_all = T)

lev0 <- unique(QTL_01$Trait)

M3 <- manhattan.plot(data = data_5, traits = lev0) + theme_classic(base_family = "Arial", base_size = 12) + theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y = element_text(size = 12), plot.tag = element_blank()) 
M3

ggsave(filename = "~/Documents/git/big_files/M0_sqrt_DS.jpg", plot = M3, width = 8, height = 8)

write.csv(QTL_02, "markers_yield_DS.csv", quote = F, row.names = F)
