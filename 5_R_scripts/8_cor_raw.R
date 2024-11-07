# Agriutilities Roza2019
rm(list = ls())
library(tidyverse)
library(metan)
library(tidyverse)
library(RColorBrewer)
library(grid)
library(gtable)
library(emmeans)
library(svglite)
library(multcomp)
library(multcompView)
library(asreml)
library(agriutilities)


############
setwd("~/Documents/git/big_files/")

a1 <- read.csv("Roza_2019_yield_raw.csv")
head(a1)
a1.2 <- a1 %>% dplyr::select("acc", colnames(a1)[9:26])
head(a1.2)

a1.1 <- a1[,c(2,4:8)] # position
head(a1.2)


lev1 <- colnames(a1)[9:26]
lev1

a2 <- a1[,c(2,9:26)] # yield

a2 <- gather(a2, key = "harv", value = "yield", 2:ncol(a2))
colnames(a2)

# 1 pound = 453.5923 g

a2$harv <- as.factor(a2$harv)
summary(a2$harv)
a2 <- a2 %>% mutate(harv = factor(harv)) %>% mutate(harv = fct_relevel(harv, lev2))
a2$g_yield <- a2$yield * 453.5923
head(a2)

a2 <- a2 %>% dplyr::select(plot, harv, g_yield)
head(a2)
a2 <- a2 %>% dplyr::filter(!harv%in% c("total_20","total_21","total_22","total_23"))

a6 <- a2 %>% spread(key = harv, value = g_yield) %>% column_to_rownames("plot")
head(a6)

cor1 <- cor(a6, use = "complete.obs")

setwd("~/Documents/git/Roza_2019/5_R_scripts/")
setwd("~/medin297@umn.edu - Google Drive/My Drive/Roza2019/cor/")
write.csv(cor1, "cor1_raw.csv", quote = F, row.names = T)

a7 <- inner_join(a1.1, a2, by = "plot")
head(a7)
str(a7)
# a7 <- a7 %>% spread(key = harv, value = sqrt_yield)

lev3 <- colnames(a7)[1:7]
a7[,lev3] <- lapply(a7[,lev3], factor)

# a8 <- na.omit(a7)
levels(a7$harv)

a9 <- a7 %>% dplyr::filter(harv %in% c("may_20", "sep_20", "may_21", "sep_21", "jun_22", "aug_22", "may_23", "jul_23"))
a9 <- droplevels(a9)

a9$water <- dplyr::recode_factor(a9$harv, 
                          "may_20" = "1" ,"sep_20" = "0",
                          "may_21" = "1" ,"sep_21" = "0",
                          "jun_22" = "1" ,"aug_22" = "0",
                          "may_23" = "1" ,"jul_23" = "0")

a9 <- a9 %>% separate(col = harv, into = c("month", "year"), sep = "_", remove = F)
str(a9)
a9$year <- as.factor(a9$year)
a9$year <- dplyr::recode_factor(a9$year, 
                                 "20" = "2020" ,"21" = "2021",
                                 "22" = "2022" ,"23" = "2023")

a10 <- na.omit(a9)

a10 <- a10[order(a10$water), ]

model2f<-asreml(fixed = HT~Test+Test:Rep,
                random = ~ at(Test):Rep:Iblock + fa(Test,1,init=initg):Genotype,
                residual = ~ dsum(~units|Test),
                data=datam)

model2f<-asreml(fixed = g_yield ~ water + water:year,
                random = ~ at(water):rep:block + fa(water,1):acc_num,
                residual = ~ dsum(~units|water),
                data=a10)
