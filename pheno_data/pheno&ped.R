rm(list = ls())

library(tidyverse)

setwd("~/Documents/git/Roza_2019/pheno_data/")
ped1 <- read.csv("Pedigree_Roza2019_1.csv")
id3 <- read.csv("Roza_ID3.csv")
id3$Plant_ID <- as.character(id3$Plant_ID)
ped1$Plant_ID <- as.character(ped1$Plant_ID)

head(id3)
head(ped1)

id4 <- id3 %>% dplyr::select(Roza2019_VCF_ID, Plant_ID)
ped2 <- ped1 %>% dplyr::select(Plant_ID, BC_ID, Female, Male)
ped3 <- inner_join(ped2, id4, by = "Plant_ID") %>% dplyr::select(-Plant_ID)

setwd("~/Documents/git/big_files/")
pheno <- read.csv("BLUE_Yi_sqrt_SpATS_DArT.csv") # ST1 (SpATS)
pheno1 <- read.csv("BLUP_Yi_sqrt_SpATS_DArT.csv")
colnames(pheno)
head(pheno1)
pheno <- pheno %>% dplyr::select(-c("total_20", "total_21", "total_22", "total_23"))
pheno1 <- pheno1 %>% dplyr::select(-c("yi_year"))

pheno2 <- inner_join(pheno, pheno1, by = c("Roza2019_VCF_ID", "PC1", "PC2", "PC3"))
colnames(pheno2)

pheno2 <- pheno2 %>% relocate(PC3, .after = yi_2023)
pheno2 <- pheno2 %>% relocate(PC2, .after = yi_2023)
pheno2 <- pheno2 %>% relocate(PC1, .after = yi_2023)
pheno2 <- pheno2 %>% relocate(yi_hrv, .after = yi_2023)

pheno3 <- inner_join(pheno2, ped3, by = "Roza2019_VCF_ID")
colnames(pheno3)

setwd("~/Documents/git/big_files/")
write.csv(pheno3, file = "BLUE&BLUP_Yi_sqrt_ped.csv", quote = F, row.names = F)


library(heatmaply)

setwd("~/Documents/git/big_files/")
# write.csv(pheno3, file = "BLUE&BLUP_Yi_sqrt_ped.csv", quote = F, row.names = F)
pheno3 <- read.csv("BLUE&BLUP_Yi_sqrt_ped.csv")

pheno4 <- pheno3[2:19]
colnames(pheno3)
G6 <- as.matrix(pheno4)
rownames(G6) <- pheno3$Roza2019_VCF_ID

?heatmaply
x <- heatmaply(G6, dendrogram = "row")
x
# y <- heatmaply(G6, cellnote = G6, cellnote_textposition = "middle center", cellnote_size = 14,)

setwd("~/medin297@umn.edu - Google Drive/My Drive/DAl21-6679/2024/10_fig/")
orca(y, "heat2.svg", width = 3, height = 3)

G6 <- as.matrix(G6)

?heatmaply
x <- heatmaply(G6)

setwd("~/medin297@umn.edu - Google Drive/My Drive/DAl21-6679/2024/10_fig/")
orca(x, "heat1.svg")