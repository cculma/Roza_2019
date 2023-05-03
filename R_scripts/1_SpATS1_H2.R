# BLUP values of yield Roza2019
rm(list = ls())

library(SpATS)
library(tidyverse)
library(data.table)
library(ggcorrplot)
library(asreml)
library(asremlPlus)
library(patchwork)

############

a1 <- read.csv("~/Documents/git/Roza_2019/raw_data/Roza2019_yield.csv")
head(a1)

colnames(a1)
lev2 <- c("acc","acc_num","rep")
a1[,lev2] <- lapply(a1[,lev2], factor)

nlevels(a1$acc)
nlevels(a1$acc_num)
a1$R <- as.factor(a1$row)
a1$C <- as.factor(a1$col)
nlevels(a1$C)
nlevels(a1$R)

lev3 <- colnames(a1)[7:23]

Y1 <- list()
for (i in 1:length(lev3)) {
  m1 <- SpATS(response = lev3[i],
              spatial = ~PSANOVA(col, row, nseg = c(nlevels(a1$C),nlevels(a1$R)), degree = c(3,3), nest.div = 2),
              fixed = ~rep,
              random = ~ rep:C + rep:R,
              genotype.as.random = T,
              genotype = "acc_num",
              data = a1,
              control = list(tolerance = 1e-03))
  
  pred4 <- as.data.frame(getHeritability(m1))
  Y1[[length(Y1)+1]] <- pred4
}

names(Y1) <- lev3


Y2 <-rbindlist(Y1, use.names=TRUE, fill=TRUE, idcol="env")
head(Y2)
colnames(Y2)[2] <- "H2"

setwd("~/Documents/git/Roza_2019/pheno_data/")
write.csv(Y2, "SpATS_H2.csv")
