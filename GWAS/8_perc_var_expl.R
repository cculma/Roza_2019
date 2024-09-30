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
trait1 <- colnames(pheno)[1:(length(colnames(pheno))-7)]
trait1

params <- set.params(fixed=c("PC1","PC2","PC3"),
                     fixed.type=rep("numeric",3), n.PC = 3)
models_1 <- c("general", "additive", "1-dom", "2-dom")

data_1 <- read.GWASpoly(ploidy=4, 
                        pheno.file="BLUE_Yi_sqrt_SpATS_DArT.csv", 
                        geno.file="Roza2019_06_GWASPoly.txt", 
                        format="numeric", n.traits=length(trait1), delim=",")


data_2 <- set.K(data = data_1, LOCO = F, n.core = 16)
data_3 <- GWASpoly(data = data_2, models = models_1, traits = trait1, params = params, n.core = 16)

# data_3 <- GWASpoly(data = data_2, models = models_1, traits = "may_20", params = params, n.core = 16)
# data_4 <- set.threshold(data_3, method= "Bonferroni", level=0.05)
# data_5 <- get.QTL(data_4)


# r^2 ---------------------------------------------------------------------
# phenotypic variance explained (PVE)

# 
# cc <- count(data_5,Trait)
# lev4 <- cc$Trait
# lev4
# cc1 <- count(data_5, Model)
# cc1$Model
# 
# data_6 <- data_5 %>% dplyr::filter(!Model %in% c("diplo-general", "diplo-additive"))
# 
# 
# data_6 <- fit.QTL(data=data_3, trait = "may_20",
#                   qtl=data_5[,c("Marker","Model")])
# 
# data_6 <- data_6 %>% group_by(Marker) %>% top_n(1, abs(R2)) %>% unite(col = "Marker1", 2:3, sep = "_", remove = T) %>% dplyr::select(Marker, Marker1, R2, pval)
# 


lev5 <- list()
for (i in 1:length(trait1)) {
  
  data_3 <- GWASpoly(data = data_2, models = models_1, traits = trait1[i], params = params, n.core = 16)
  data_4 <- set.threshold(data_3, method= "Bonferroni", level=0.05)
  data_5 <- get.QTL(data_4)
  if (nrow(data_5) == 0) next
  data_6 <- fit.QTL(data=data_3, trait = trait1[i], qtl=data_5[,c("Marker","Model")])
  
  data_7 <- data_6 %>% group_by(Marker) %>% top_n(1, abs(R2)) %>% unite(col = "Marker1", 2:3, sep = "_", remove = T) %>% dplyr::select(Marker, Marker1, R2, pval)
  
  data_7$trait <- trait1[i]
  lev5[[length(lev5) + 1]] <- data_7
}
