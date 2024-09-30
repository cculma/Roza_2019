# ID match
library(tidyverse)

setwd("~/Documents/git/big_files/")

PCA <- read.csv("PCA_Roza2019.csv")
PCA <- PCA %>% dplyr::select(Roza2019_VCF_ID, Sample_ID)
colnames(PCA)

GEBV1 <- read.csv("GEBV_GBLUP_st1.csv")
GEBV2 <- read.csv("GEBV_GBLUP_st2.csv")
colnames(GEBV1)[1] <- "Roza2019_VCF_ID"
colnames(GEBV2)[1] <- "Roza2019_VCF_ID"

GEBV1 <- inner_join(PCA, GEBV1, by = "Roza2019_VCF_ID") %>% dplyr::select(-Roza2019_VCF_ID)
GEBV2 <- inner_join(PCA, GEBV2, by = "Roza2019_VCF_ID") %>% dplyr::select(-Roza2019_VCF_ID)

write.csv(GEBV1, "GEBV_GBLUP_st1.1.csv", quote = F, row.names = F)
write.csv(GEBV2, "GEBV_GBLUP_st2.1.csv", quote = F, row.names = F)


BLUES<- read.csv("BLUE_Yi_sqrt_SpATS_DArT.csv")
BLUPS <- read.csv("BLUP_Yi_sqrt_SpATS_DArT.csv")

BLUES <- inner_join(PCA, BLUES, by = "Roza2019_VCF_ID") %>% dplyr::select(-Roza2019_VCF_ID)
BLUPS <- inner_join(PCA, BLUPS, by = "Roza2019_VCF_ID") %>% dplyr::select(-Roza2019_VCF_ID)


write.csv(BLUES, "BLUE_Yi_sqrt_SpATS_DArT.1.csv", quote = F, row.names = F)
write.csv(BLUPS, "BLUP_Yi_sqrt_SpATS_DArT.1.csv", quote = F, row.names = F)
