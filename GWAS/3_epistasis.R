# MAEIF1 prepare matrix

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
library(GWASpoly)

setwd("~/Documents/git/big_files/")

load("~/Documents/git/big_files/Yi_st1_51081_sqrt.RData")
data_3 <- set.threshold(Yi_data_3, method= "Bonferroni", level=0.05)
QTL_03 <- get.QTL(data_3)

load("~/Documents/git/big_files/Yi_st2_51081_sqrt.RData")
data_4 <- set.threshold(Yi_data_3, method= "Bonferroni", level=0.05)
QTL_04 <- get.QTL(data_4)

load("~/Documents/git/big_files/data_Yi_DS_20.RData")
data_5 <- set.threshold(Yi_data_3, method= "Bonferroni", level=0.05)
QTL_05 <- get.QTL(data_5)

QTL_06 <- rbind(QTL_03, QTL_04, QTL_05)
QTL_06.1 <- QTL_06 %>% distinct(Marker, .keep_all = T)

QTL_06.1 <- QTL_06.1 %>% unite(col = "marker1", c("Chrom", "Position"), sep = "_", remove = F)


setwd("~/Documents/git/big_files/")
Y <- read.csv("BLUP_Yi_sqrt_SpATS_DArT.csv", header = T)
colnames(Y)
head(Y)
Y <- Y %>% dplyr::select(Roza2019_VCF_ID, yi_hrv)
Y1 <- na.omit(Y) 
Y1$yi_hrv <- (Y1$yi_hrv * 453.592)

G <- read.table("Roza2019_06_GS.txt", header = TRUE, check.names = F)
G[1:5,1:5]
dim(G) # 424 51081
common <- intersect(Y1$Roza2019_VCF_ID,rownames(G))

marks <- G[common,]
marks[1:5,1:5]
dim(marks) # 424 51081
class(marks)

marks.1 <- marks %>% dplyr::select(QTL_06.1$marker1)
dim(marks.1) # 424  137

Y2 <- Y1[match(common, Y1$Roza2019_VCF_ID),]
dim(Y2) # 424  2
Y2 <- Y2 %>% remove_rownames() %>% column_to_rownames(var = "Roza2019_VCF_ID")
colnames(Y2)
Y3 <- Y2

data <- merge(as.data.frame(Y3), as.data.frame(marks.1), by = 'row.names', all = TRUE) %>% column_to_rownames(var = 'Row.names')
data[1:5,1:5]
colnames(data)[1] <- "PHENOTYPE"
data <- na.omit(data)
dim(data) # 424  138 : 424 genotypes and 137 markers

setwd("~/Documents/git/big_files/")
write.csv(data, "epistasis01.csv", quote = F, row.names = T)


# Longitudinal Dataset ----------------------------------------------------


pheno <- read.csv("BLUE_Yi_sqrt_DArT_DS1.csv")
head(pheno)
pheno1 <- pheno %>% dplyr::select(Roza2019_VCF_ID, X23, Stress)
pheno1$Stress <- as.factor(pheno1$Stress)
levels(pheno1$Stress)
pheno1$Stress <- recode_factor(pheno1$Stress, "DS" = "0", "NS" = "1")
head(pheno1)
pheno1$X23 <- (pheno1$X23 * 453.592)
hist(pheno1$X23)

load("~/Documents/git/big_files/data_Yi_DS_20.RData")
data_5 <- set.threshold(Yi_data_3, method= "Bonferroni", level=0.05)
QTL_05 <- get.QTL(data_5)
QTL_05 <- QTL_05 %>% dplyr::filter(Trait == "X23")
QTL_05 <- QTL_05 %>% distinct(Marker, .keep_all = T)
QTL_05 <- QTL_05 %>% unite(col = "marker1", c("Chrom", "Position"), sep = "_", remove = F)


G <- read.table("Roza2019_06_GS.txt", header = TRUE, check.names = F)
G[1:5,1:5]
dim(G) # 424 51081
marks.1 <- G %>% dplyr::select(QTL_05$marker1)
dim(marks.1) # 424  40
marks.1 <- marks.1 %>% rownames_to_column("Roza2019_VCF_ID")
marks.1[1:5,1:5]

data <- left_join(pheno1, marks.1, by = "Roza2019_VCF_ID") %>% dplyr::select(-Roza2019_VCF_ID)

data[1:5,1:5]
colnames(data)[1] <- "PHENOTYPE"
data <- na.omit(data)
dim(data) # 848  42
print(summary(data[,1:15]))

setwd("~/Documents/git/big_files/")
write.csv(data, "epistasis02.csv", quote = F, row.names = T)
