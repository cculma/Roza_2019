# AMatrix

library(AGHmatrix)
library(ASRgenomics)
library(tidyverse)

setwd("~/Documents/git/Roza_2019/pheno_data/")

# rows
# 1-60	parents
# 61-420	F1
# 1-420	parents and F1
# 421-844	BC


ped1 <- read.csv("Pedigree_Roza2019.csv", row.names = 1)
head(ped1)
ped3 <- ped1[1:420,]
ped2 <- ped1[421:844,]
PCA1
ped2$Plant_ID <- as.integer(ped2$Plant_ID)
ped2 <- inner_join(PCA1, ped2, by = "Plant_ID") %>% dplyr::select(-Plant_ID)
colnames(ped2)[1] <- "Plant_ID"
head(ped2)
ped4 <- rbind(ped3, ped2)

APED1 <- Amatrix(ped4, ploidy = 4)

APED1[1:5,1:5]
ASRgenomics::kinship.heatmap(APED1)
dim(APED1)
kin1 <- ASRgenomics::kinship.heatmap(APED1[421:844,421:844]) # only BC
class(kin1)
# save pdf 8 by 8

APED2 <- as.matrix(APED1[421:844,421:844])

APED2 <- as.data.frame(APED1[421:844,421:844])
APED2[1:5,1:5]
APED2 <- APED2 %>% rownames_to_column("row")
setwd("~/Documents/git/big_files/")
write.table(APED2, "Roza2019_AMatrix.tsv", row.names = F, quote = F, sep = "\t")
?kinship.heatmap
?qc.filtering


?Amatrix
APED3 <- Amatrix(ped4, ploidy=4, slater = TRUE)
APED3 <- Amatrix(ped4, ploidy=4, w=0.1, slater = TRUE)
APED3 <- Amatrix(ped4, ploidy=4, w=0.2, slater = TRUE)
APED3 <- Amatrix(ped4, ploidy=4, w=0.3, slater = TRUE)
ASRgenomics::kinship.heatmap(APED3[421:844,421:844])
APED4 <- as.data.frame(APED3[421:844,421:844])
APED4[1:5,1:5]

# read Pheno --------------------------------------------------------------

setwd("~/Documents/git/big_files/")

PCA <- read.csv("PCA_Roza2019.csv")
colnames(PCA)
PCA1 <- PCA %>% dplyr::select(Roza2019_VCF_ID, Plant_ID)

pheno <- read.csv("BLUE_Yi_sqrt_SpATS_DArT.csv")
colnames(pheno)
pheno <- pheno[,1, drop = F]
PCA2 <- inner_join(pheno, PCA1, by = "Roza2019_VCF_ID")

# read Geno ---------------------------------------------------------------

setwd("~/Documents/git/big_files/")
GS1 <- read.table("Roza2019_06_GS.txt", row.names = 1)
# GS1 <- GS1 %>% rownames_to_column("Roza2019_VCF_ID")

GS1[1:5,1:5]
dim(GS1) # 424 51081
# GS1 <- qc.filtering(M = GS1, base = F, ref = NULL, maf = 0.05, marker.callrate = 0.2, ind.callrate = 0.2, heterozygosity = 0.9, Fis = 0.5, impute = F, plots = T)


# GS2 <- inner_join(PCA2, GS1, by = "Roza2019_VCF_ID")
# GS2[1:5,1:5]
# GS2 <- GS2 %>% dplyr::select(-Roza2019_VCF_ID) %>% column_to_rownames("Plant_ID")
# dim(GS2) # 424 51081
# GS2[1:5,1:5]
# rownames(GS2)
# class(GS2)
# GS2 <- as.matrix(GS2)

GS1 <- as.matrix(GS1)
G_Van <- Gmatrix(GS1, method="VanRaden", ploidy=4) # additive relationship matrix
G_Ful <- Gmatrix(GS1, method="Slater", ploidy=4)   # full-autopolyploid matrix
G_Dom <- Gmatrix(GS1, method="Endelman", ploidy=4) # dominance (digenic) matrix
G_Van[1:5,1:5]
G_Ful[1:5,1:5]
G_Dom[1:5,1:5]

G_Van1 <- kinship.diagnostics(G_Van)
kinship.heatmap(G_Van)
kinship.heatmap(G_Ful)
kinship.heatmap(G_Dom)

G2A_Van <- match.G2A(A = APED2, G = G_Van, clean = T, ord = T, mism = T, RMdiff = T)
G2A_Ful <- match.G2A(A = APED2, G = G_Ful, clean = T, ord = T, mism = T, RMdiff = T)
G2A_Dom <- match.G2A(A = APED2, G = G_Dom, clean = T, ord = T, mism = T, RMdiff = T)

G2A_Van$plotG2A
G2A_Ful$plotG2A
G2A_Dom$plotG2A

G_blend <- G.tuneup(G = G_Van, blend = T, pblend = 0.05)$Gb
G_bend <- G.tuneup(G = G_Van, bend = T)$Gb

kinship.heatmap(G_Van)
kinship.heatmap(G_blend)
kinship.heatmap(G_bend)

G2A_bend <- match.G2A(A = APED2, G = G_bend, clean = T, ord = T, mism = T, RMdiff = T)
G2A_bend$plotG2A


#Computing H matrix (Martini)
Hmat_Ma <- Hmatrix(A=APED2, G=G_bend, method="Martini", 
                        ploidy=4, missingValue=-9, maf=0.05)

#Computing H matrix (Munoz) 2014
Hmat_Mu <- Hmatrix(A=APED2, G=G_bend, markers = GS1, 
                      ploidy=4, method="Munoz", 
                      missingValue=-9, maf=0.05)


kinship.heatmap(Hmat_Ma)
kinship.heatmap(Hmat_Mu)

?H.inverse
