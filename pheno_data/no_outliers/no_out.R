rm(list = ls())

library(SpATS)
library(tidyverse)

setwd("~/Documents/git/Roza_2019/pheno_data/no_outliers/")

a1 <- read.csv("ST0_sep_21.csv")
a2 <- read.csv("ST0_aug_22.csv")
a1 <- a1 %>% dplyr::select(1,2)
a2 <- a2 %>% dplyr::select(1,2)

a3 <- inner_join(a1,a2, by = 'acc')
head(a3)

min(a3$ST0_sep_21)
max(a3$ST0_sep_21)
mean(a3$ST0_sep_21)

min(a3$ST0_aug_22)
max(a3$ST0_aug_22)
mean(a3$ST0_aug_22)

setwd("~/Documents/git/Roza_2019/pheno_data/")
Y1 <- read.csv("deglos.csv")
Y1 <- Y1[,c(1:2)]
colnames(Y1) <- c("acc", "gen")
Y1$gen <- as.character(Y1$gen)

setwd("~/Documents/git/Roza_2019/pheno_data/")
Y2 <- read.csv("VCF_IDs_Roza2019.csv", header = T, check.names = F)
Y2 <- Y2[,c(2,4)]


setwd("~/Documents/git/Roza_2019/pheno_data/")
Y3 <- read.csv("PCA_Roza2019.csv")

Y4 <- inner_join(Y1, Y2, by = "gen") %>% inner_join(., a3, by = "acc")%>% inner_join(., Y3, by = "VCF")
Y4 <- Y4[,-c(1,2)]

setwd("~/Documents/git/Roza_2019/pheno_data/no_outliers/")
write.csv(Y4, "no_out1.csv", quote = F, row.names = F)



setwd("~/Documents/git/Roza_2019/pheno_data/")

a3 <- read.csv("BLUPs_Yi_Roza2019.csv")

min(a3$ST0_sep_21)
max(a3$ST0_sep_21)
mean(a3$ST0_sep_21)

min(a3$ST0_aug_22)
max(a3$ST0_aug_22)
mean(a3$ST0_aug_22)
