rm(list = ls())
library(tidyr)
library(tidyverse)
library(sommer)
library(biganalytics)
library(bigmemory)
library(parallel)
library(doParallel)
library(iterators)
library(foreach)
library(gbm)
library(rrBLUP)
library(e1071)
library(caret)
library(dplyr)
library(ranger)
library(nnet)
library(brnn)
library(arm)
library(monomvn)
library(GROAN)
library(AGHmatrix)
library(ggfortify)
library(ggplot2)
library(ASRgenomics)
library(factoextra)
library(FactoMineR)

cl <- makePSOCKcluster(15)
registerDoParallel(cl)

load("~/OneDrive - Washington State University (email.wsu.edu)/Roza_2019/GS1/PCA_Roza2019.R")
load("/Users/cesarmedina/Library/CloudStorage/OneDrive-WashingtonStateUniversity(email.wsu.edu)/Roza_2019/GS1/PCA_Roza2019.RData")
# load("/home/hawkins/Documents/Cesar/blup_data/Roza2019/git/r_data/PCA_Roza2019.RData")
# load("/home/hawkins/Documents/Cesar/blup_data/Roza2019/git/r_data/PCA_Roza2019.RData")


###########################
# geno data

G <- read.csv("~/Documents/git/big_files/Roza2019_05_GWASPoly.txt", header = TRUE, row.names = 1, check.names = F)
# G <- read.csv("~/Documents/blup_data/Roza2019/Analysis_2021/GWAS/MPP_Ms2_GWASPoly.txt", header = TRUE, row.names = 1, check.names = F)
G[1:5,1:5]
G1 <- G %>% unite(Chrom1, 1:2, remove = T)
G1 <- as.matrix(G1 %>% remove_rownames() %>% column_to_rownames(var = "Chrom1"))
G2 <- t(G1)
G2[1:5,1:5]
G2 <- as.data.frame(G2)
str(G2)
numo <- atcg1234(data=G2, ploidy=4); 
numo$M[1:5,1:5]
numo$ref.allele[,1:5]

head(Y5$Roza2019_VCF_ID)
have.both = intersect(rownames(numo$M), Y5$Roza2019_VCF_ID)
length(have.both)

G3 <- numo$M[have.both,]
G3[1:5,1:5]
dim(G3) # 424 62839
class(G3)
geno.ps <- G3 / 4
geno.ps[1:5,1:5]
Mpr <- snp.pruning(M = geno.ps, pruning.thr = 0.90, window.n = 100, by.chrom = F, overlap.n = 10, seed = 1208, iterations = 20)
G4 <- Mpr$Mpruned * 4
G4[1:5,1:5]
dim(G4) # 424 51081
G5 <- as.data.frame(t(G4))
G5[1:5,1:5]
G5 <- G5 %>% tibble::rownames_to_column("snp_name") %>% separate(col = 1, into = c("Chrom", "Position"), remove = T, sep = "_")
Marker <- seq(1:nrow(G5))
G6 <- cbind(Marker, G5)
G6[1:5,1:5]

setwd("~/Documents/git/big_files/")
write.csv(G6, "Roza2019_06_GWASPoly.txt", row.names = F, quote = F)

# NZV <- nearZeroVar(G4)
# G5 <- G4[,NZV]
MAF <- as.data.frame((colSums(G4)/ 424) / 4) %>% rownames_to_column("marker")


G_OK <- Gmatrix(G4, method = "Slater", ploidy=4)
det(G_OK)
G_OK[1:5,1:5]
G_OK_Inv <- solve(G_OK)



Y1 <- Y5[match(have.both,Y5$Roza2019_VCF_ID),]
dim(Y1)

PCs <- prcomp(G4, scale=F, center = T)
PC2s <- as.data.frame(PCs$x[,1:3])
autoplot(PCs) + theme_bw()
fviz_eig(PCs, addLabels = T, ylim = c(0, 4))


PC1s <- prcomp(G_OK, scale=F, center = T)
autoplot(PC1s) + theme_bw()
fviz_eig(PC1s, addLabels = T, ylim = c(0, 4))

PC2s <- as.data.frame(PCs$x[,1:3]) %>% rownames_to_column("Roza2019_VCF_ID")

######################
# pheno
# Y <- read.csv("~/Documents/git/Roza_2019/pheno_data/pheno_nph.csv", header = T)
Y <- read.csv("~/Documents/git/Roza_2019/pheno_data/BLUE_Yi_Roza2019.csv", header = T)
colnames(Y)[1] <- "Plant_ID"
Y$Plant_ID <- as.character(Y$Plant_ID)

Y1 <- read.csv("~/Documents/git/Roza_2019/pheno_data/Roza_ID3.csv", header = T)
head(Y1)

Y5 <- inner_join(Y1, Y, by = "Plant_ID")
head(Y5)
dim(Y5)

## https://cran.r-project.org/web/packages/ggfortify/vignettes/plot_pca.html

Y6 <- inner_join(Y5, PC2s, by = "Roza2019_VCF_ID")
head(Y6)
Y7 <- Y6 %>% dplyr::select(-c("Roza2019_VCF_ID1","Sample_ID","Plant_ID"))
setwd("~/Documents/git/big_files/")
write.csv(Y7, "Yi_st1.csv", row.names = F, quote = F)


Y <- read.csv("~/Documents/git/Roza_2019/pheno_data/BLUP_Yi_Roza2019.csv", header = T)
colnames(Y)[1] <- "Plant_ID"
Y$Plant_ID <- as.character(Y$Plant_ID)
Y1 <- read.csv("~/Documents/git/Roza_2019/pheno_data/Roza_ID3.csv", header = T)
head(Y1)

Y5 <- inner_join(Y1, Y, by = "Plant_ID")
head(Y5)
dim(Y5)
Y6 <- inner_join(Y5, PC2s, by = "Roza2019_VCF_ID")
head(Y6)
Y7 <- Y6 %>% dplyr::select(-c("Roza2019_VCF_ID1","Sample_ID","Plant_ID"))
head(Y7)

setwd("~/Documents/git/big_files/")
write.csv(Y7, "Yi_st2.csv", row.names = F, quote = F)


classif <- kmeans(x = PC2s[,1:2], centers = 3, iter.max = 1000, nstart = 10, algorithm = "Lloyd", trace = F)
PC2s$group <- classif$cluster

p1 <- ggplot(PC2s, aes(x = PC1, y = PC2, color = group)) + geom_point() + theme_minimal() 
p1

nb.cols <- 15
mycolors <- colorRampPalette(brewer.pal(8, "Paired"))(nb.cols)
