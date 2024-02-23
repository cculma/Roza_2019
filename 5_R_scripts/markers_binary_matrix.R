# this script allows to genrate a binary table which summarize redundant markers in non-redundant markers in rows and GWASpoly model or trait in columns.

rm(list = ls()) # clean Global Environment
# setwd("~/Documents/Cesar/blup_data/Roza2019/Analysis_2021/GWAS/")
# setwd("~/OneDrive - Washington State University (email.wsu.edu)/Roza_2019/git/Roza2019/")
library(GWASpoly)
library(tidyverse)
library(reshape2)

load("~/OneDrive - Washington State University (email.wsu.edu)/Roza_2019/GWASpoly_results/data_3_80177_year.RData")


S1 <- QTL_01 %>% dplyr::select(2,4) %>% distinct(Marker, Model, .keep_all = T) 
S1 <- dcast(S1, formula = Marker ~ Model, fun.aggregate = length)

S2 <- QTL_01 %>% dplyr::select(1,4) %>% distinct(Marker, Trait, .keep_all = T) 
S2 <- dcast(S2, formula = Marker ~ Trait, fun.aggregate = length)
colnames(S2)
S4 <- S2[,c(1,11,9,7,2,22,12,10,8,3,23,18,17,16,19,13,14,5,15,20,21,6,4)]
S3 <- inner_join(QTL_06, S1, by = "Marker") %>% inner_join(., S4, by = "Marker")

write.table(S3, "~/OneDrive - Washington State University (email.wsu.edu)/Roza_2019/GWASpoly_results/S3.tsv", row.names = F, quote = F, sep = "\t")


# to generate frequency of markers shared by trait

S2 <- QTL_01 %>% filter(Trait != "MET_all") %>% dplyr::select(1,4) %>% distinct(Marker, Trait, .keep_all = T) 

cc <- count(S2, Marker)
cc1 <- count(cc, n)

fig <- plot_ly(
  x = cc1$n,
  y = cc1$nn,
  name = "freq Markers",
  type = "bar", text = cc1$nn, textposition = 'auto') %>% layout(xaxis = list(autotypenumbers = 'strict', title = 'Trait'), yaxis = list(title = 'Freq')) %>% config(toImageButtonOptions = list(format = "svg",filename = "fig", width = 600, height = 300))

orca(fig, "bar1.svg", width = 5 * 96, height = 4 * 96)
