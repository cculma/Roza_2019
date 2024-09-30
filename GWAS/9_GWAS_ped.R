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

pheno <- read.csv("BLUE&BLUP_Yi_sqrt_ped.csv", row.names = 1) # BLUE, BLUP & PED 

head(pheno)
trait1 <- colnames(pheno)[1:(length(colnames(pheno))-6)]
trait1
params <- set.params(fixed=c("PC1","PC2","PC3","BC_ID","Female","Male"),
                     fixed.type= c("numeric","numeric","numeric","factor","factor","factor"), n.PC = 3)
models_1 <- c("general", "additive", "1-dom", "2-dom",  "diplo-additive", "diplo-general")

# numeric
# Roza2019_06_GWASPoly.txt # DArTag

data_1 <- read.GWASpoly(ploidy=4, 
                        pheno.file="BLUE&BLUP_Yi_sqrt_ped.csv", 
                        geno.file="Roza2019_06_GWASPoly.txt", 
                        format="numeric", n.traits=length(trait1), delim=",")

data_2 <- set.K(data = data_1, LOCO = T, n.core = 16)

Ped_data_3 <- GWASpoly(data = data_2, models = models_1, traits = "may_23", params = params, n.core = 16)

data_4 <- set.threshold(Ped_data_3, method= "Bonferroni", level=0.05)
QTL_04 <- get.QTL(data_4)


manhattan.plot(data = data_4, traits = "may_23") + theme_classic(base_family = "Arial", base_size = 12) + theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y = element_text(size = 12), plot.tag = element_blank()) 

# ggsave(filename = "~/Documents/git/big_files/M0_sqrt.jpg", plot = M0, width = 16, height = 16)

# compare with previous results  ------------------------------------------

# load("~/Documents/git/big_files/Yi_st1_51081_sqrt.RData")
data_3 <- set.threshold(Yi_data_3, method= "Bonferroni", level=0.05)
QTL_03 <- get.QTL(data_3)

cc <- dplyr::count(QTL_03, Trait)
QTL_03.1 <- QTL_03 %>% dplyr::filter(Trait == "may_23")

manhattan.plot(data = data_3, traits = "may_23") + theme_classic(base_family = "Arial", base_size = 12) + theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y = element_text(size = 12), plot.tag = element_blank()) 


# r^2 ---------------------------------------------------------------------
# phenotypic variance explained (PVE)R
data_2 <- set.K(data = data_1, LOCO = F, n.core = 16)

Ped_data_3.1 <- GWASpoly(data = data_2, models = models_1, traits = "may_23", params = params, n.core = 16)


data_4 <- set.threshold(Ped_data_3.1, method= "Bonferroni", level=0.05)
QTL_04 <- get.QTL(data_4)

manhattan.plot(data = data_4, traits = "may_23") + theme_classic(base_family = "Arial", base_size = 12) + theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y = element_text(size = 12), plot.tag = element_blank()) 

if (nrow(data_5) == 0) next


QTL_05 <- QTL_04 %>% dplyr::filter(!Model %in% c("diplo-general", "diplo-additive"))

data_6 <- fit.QTL(data=data_4, trait = "may_23", qtl=QTL_05[,c("Marker","Model")])

data_7 <- data_6 %>% group_by(Marker) %>% top_n(1, abs(R2)) %>% unite(col = "Marker1", 2:3, sep = "_", remove = T) %>% dplyr::select(Marker, Marker1, R2, pval)



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
