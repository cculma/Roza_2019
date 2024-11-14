library(tidyverse)

setwd("~/Documents/git/big_files/")

PCA <- read.csv("PCA_Roza2019.csv")
head(PCA)

PCA1 <- PCA %>% dplyr::select( Susceptible_parent1, Roza2019_VCF_ID)
head(PCA1)
PCA1$Susceptible_parent1 <- paste0("P",PCA1$Susceptible_parent1)


 
# PCA$pop2 <- recode_factor(PCA$Susceptible_parent1, 
#                           "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8","P9",
#                           "P10", "P11", "P12", "P13", "P14", "P15", "P16","P17",
#                           "P18", "P19", "P20", "P21", "P22", "P23", "P24", "P25", "P26", "P27", "P28",
#                           "P29", "P30")
                          
cc06 <- dplyr::count(PCA1, Susceptible_parent1)
cc06$Susceptible_parent1
cc06$seq <- seq(1:30)
cc07 <- inner_join(cc06, PCA1, by = "Susceptible_parent1")
cc07 <- cc07 %>% dplyr::select(seq, Roza2019_VCF_ID)


G1 <- read.table("Roza2019_06_GS.txt", sep = "\t", row.names = 1)
G1[1:5,1:5]
str(G1)
dim(G1)
G1 <- G1 %>% mutate_all(as.character)
G1 <- as.matrix(G1)

G1[G1 == "0"] <- "001001001001"
G1[G1 == "1"] <- "001001001002"
G1[G1 == "2"] <- "001001002002"
G1[G1 == "3"] <- "001002002002"
G1[G1 == "4"] <- "002002002002"

df4 <- as.data.frame(G1)
df4 <- df4 %>% rownames_to_column("Roza2019_VCF_ID")
df4[1:5,1:5]
dim(df4)

col9.2 <- inner_join(cc07, df4, by = "Roza2019_VCF_ID")
colnames(col9.2)[1] <- "pop"
colnames(col9.2)[2] <- "ind"
col9.2[1:10, 1:5]

setwd("~/Documents/git/big_files/")
write.table(col9.2, "Roza2019_06_genodive.txt", sep = "\t", col.names = T, row.names = F, quote = F) # 2247 2999
