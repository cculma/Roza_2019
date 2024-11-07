# env value

rm(list = ls())
library(tidyverse)
library(ggpubr)
library(svglite)

setwd("~/Documents/git/big_files/")
h2 <- read.csv("H2_h2g.csv")
pheno <- read.csv("BLUE_Yi_sqrt_SpATS_DArT.csv")
colnames(pheno)
pheno <- pheno[,1:19]
pheno <- pheno %>% column_to_rownames("Roza2019_VCF_ID")

pheno1 <- pheno %>% rownames_to_column("Roza2019_VCF_ID")
pheno1 <- pheno1 %>% dplyr::filter(Roza2019_VCF_ID %in% c("423_S314_R1_001", "271_S148_R1_001", "584_S363_R1_001", "513_S322_R1_001", "24_S21_R1_001"))
pheno1 <- pheno1 %>% column_to_rownames("Roza2019_VCF_ID")
pheno1 <- as.data.frame(t(pheno1))
pheno1 <- pheno1 %>% rownames_to_column("env")

env1 <- as.data.frame(colMeans(pheno))
colnames(env1)[1] <- "env_val"
env1 <- env1 %>% rownames_to_column("env")

env2 <- inner_join(env1, h2, by = "env") %>% inner_join(., pheno1, by = "env")
head(env2)
env2$env_val <- env2$env_val * 453.592

env3 <- inner_join(env1, pheno1, by = "env")
env3 <- env3 %>% gather(key = "Roza2019_VCF_ID", value = "BLUE", 3:ncol(env3))
env3$BLUE <- env3$BLUE * 453.592


ggscatter(env2, x = "env_val", y = "H2", add = "reg.line", fullrange = T, rug = T)

plot1 <- ggscatter(env3, x = "env_val", y = "BLUE", add = "reg.line", color = "Roza2019_VCF_ID", palette = "jco", alpha = 0.4, fullrange = F, rug = T, point = T,) + stat_cor(aes(color = Roza2019_VCF_ID, label = after_stat(rr.label)), label.x = 0.4) + stat_regline_equation(aes(color = Roza2019_VCF_ID), label.x = 0.5) + stat_cor(aes(color = Roza2019_VCF_ID)) 

plot1 <- ggscatter(env3, x = "env_val", y = "BLUE", add = "reg.line", color = "Roza2019_VCF_ID", palette = "jco", alpha = 0.4, fullrange = F, rug = T, point = T,) + stat_cor(aes(color = Roza2019_VCF_ID, label = after_stat(rr.label)), label.x = 0.4) + stat_regline_equation(aes(color = Roza2019_VCF_ID), label.x = 0.5)


setwd("~/medin297@umn.edu - Google Drive/My Drive/Roza2019/figs/")
svglite(filename = "env_value.svg", width = 4.5, height = 4, fix_text_size = F)
plot(plot1)
invisible(dev.off())
