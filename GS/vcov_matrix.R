setwd("~/Documents/git/Roza_2019/GS/")
rm(list = ls())

library(tidyverse)
library(data.table)
library(reshape2)

setwd("~/Documents/git/big_files/")
BLUE <- read.csv("BLUE_Yi_sqrt_SpATS_DArT.csv", row.names = 1)
BLUP <- read.csv("BLUP_Yi_sqrt_SpATS_DArT.csv", row.names = 1)

head(BLUE)
BLUE <- BLUE %>% dplyr::select(-c(PC1, PC2, PC3))
BLUP <- BLUP %>% dplyr::select(-c(PC1, PC2, PC3))


BLUE <- BLUE %>% dplyr::select(-c(total_20, total_21, total_22, total_23))


# colnames(BLUE)
# class(BLUE)
# BLUE <- as.matrix(BLUE)
a1 <- cov(BLUE)
a2 <- cov(BLUP)

a1[upper.tri(a1)] <- NA
# diag(a1) <- NA
# a1 <- reshape2::melt(a1, na.rm = T)
# a1$value <- round(a1$value, 3)

a2[upper.tri(a2)] <- NA
# diag(a2) <- NA
# a2 <- reshape2::melt(a2, na.rm = T)
# a2$value <- round(a2$value, 3)

colnames(a2)
a4 <- as.data.frame(diag(a2))

# rowMeans 
# a2[upper.tri(a2)] <- NA
diag(a2) <- NA
a3 <- as.data.frame(rowMeans(a2, na.rm = T))

a3$`rowMeans(a2, na.rm = T)` <- round(a3$`rowMeans(a2, na.rm = T)`, 3)

setwd("~/medin297@umn.edu - Google Drive/My Drive/Roza2019/cor/")
b1 <- read.csv("cor_st1_sqrt.csv", row.names = 1)
b2 <- read.csv("cor_st2_sqrt.csv", row.names = 1)


cor1 <- matrix(NA, nrow = nrow(b1), ncol = ncol(b1))
colnames(cor1) <- rownames(cor1) <- colnames(b1)
cor1[upper.tri(cor1)] <- a1[upper.tri(a1)]
diag(cor1) <- diag(a1)
cor1[lower.tri(cor1)] <- b1[lower.tri(b1)]


cor2 <- matrix(NA, nrow = nrow(b2), ncol = ncol(b2))
colnames(cor2) <- rownames(cor2) <- colnames(b2)
cor2[upper.tri(cor2)] <- a2[upper.tri(a2)]
diag(cor2) <- diag(a2)
cor2[lower.tri(cor2)] <- b2[lower.tri(b2)]

setwd("~/medin297@umn.edu - Google Drive/My Drive/Roza2019/cor/")
write.csv(cor1, "vcov_ST1.csv", quote = F, row.names = T)
write.csv(cor2, "vcov_ST2.csv", quote = F, row.names = T)


setwd("~/medin297@umn.edu - Google Drive/My Drive/Roza2019/cor/")
cor1 <- read.csv("cor_st1_sqrt.csv", row.names = 1)
colnames(cor1)
cor1 <- cor1[1:18,1:18]


cor1[upper.tri(cor1)] <- NA
diag(cor1) <- NA

cor1 <- cor1 %>% rownames_to_column("col1")
cor1 <- melt(cor1, na.rm = T, id.vars = 1)
melt(d, id.vars="a")
cor1$value <- round(cor1$value, 3)

