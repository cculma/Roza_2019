
rm(list = ls())

library(heatmaply)
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


head(Yi_data_3@map)
map1 <- Yi_data_3@map
map1 <- map1 %>% unite(col = "marker", c("Chrom","Position"), sep = "_")

ggscatter(score3, x = "may_20", y = "jun_20", conf.int = T, point = T, add = "reg.line") + stat_regline_equation(label.x = 3, label.y = 8)


lev1 <- names(Yi_data_3@scores)
lev1 <- lev1[1:18]

score2 <- data.frame(matrix(nrow = 51081))

for (i in 1:length(lev1)) {
  
  score1 <- Yi_data_3@scores[[lev1[i]]]
  score1 <- score1 %>% dplyr::select(general)
  rownames(score1) <- map1$marker
  colnames(score1)[1] <- lev1[i]
  score2 <- cbind(score2, score1)
  
}
score2 <- score2[,-1]

score3 <- score2

score3[score3 < 0.1 ] = NA

G6 <- cor(score3, use = "complete.obs")

diag(G6) <- NA

x <- heatmaply(G6)
x


setwd("~/medin297@umn.edu - Google Drive/My Drive/DAl21-6679/2024/10_fig/")
orca(y, "heat2.svg", width = 3, height = 3)

G6 <- as.matrix(G6)
