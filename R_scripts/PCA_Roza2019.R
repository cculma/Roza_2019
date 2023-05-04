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

# install.packages(c("GROAN","AGHmatrix","ggfortify"))

cl <- makePSOCKcluster(15)
registerDoParallel(cl)

load("~/OneDrive - Washington State University (email.wsu.edu)/Roza_2019/GS1/PCA_Roza2019.R")
load("/Users/cesarmedina/Library/CloudStorage/OneDrive-WashingtonStateUniversity(email.wsu.edu)/Roza_2019/GS1/PCA_Roza2019.RData")
# load("/home/hawkins/Documents/Cesar/blup_data/Roza2019/git/r_data/PCA_Roza2019.RData")
# load("/home/hawkins/Documents/Cesar/blup_data/Roza2019/git/r_data/PCA_Roza2019.RData")
######################
# pheno
Y <- read.csv("~/Documents/Cesar/blup_data/Roza2019/git/Roza2019/pheno_data/pheno_nph.csv", header = T)
Y1 <- read.csv("~/Documents/Cesar/blup_data/Roza2019/git/Roza2019/pheno_data/deglos.csv", header = T)
colnames(Y1)[2] <- "Acc_roza2019"
Y1$Acc_roza2019 <- as.factor(Y1$Acc_roza2019)
head(Y1)
Y2 <- read.csv("/home/hawkins/Documents/Cesar/blup_data/Roza2019/git/Roza2019/pheno_data/VCF_IDs.csv", header = T, check.names = F)
Y2 <- Y2[,c(1,3,8)]
Y3 <- inner_join(Y1, Y2, by = "Acc_roza2019")
head(Y3)
Y3 <- Y3[,-4]
factor(Y3$Male_parent_number)
factor(Y3$Acc_roza2019)
cc <- count(Y3, Planting_ID)
## https://cran.r-project.org/web/packages/ggfortify/vignettes/plot_pca.html
Y4 <- Y3  %>% column_to_rownames(var = 'VCF')

###########################
# geno data

G <- read.csv("~/Documents/Cesar/blup_data/Roza2019/Analysis_2021/GWAS/MPP_Ms2_GWASPoly.txt", header = TRUE, row.names = 1, check.names = F)

G1 <- G %>% unite(Chrom1, 1:2, remove = T)
G1 <- as.matrix(G1 %>% remove_rownames() %>% column_to_rownames(var = "Chrom1"))
G2 <- t(G1)
G2[1:5,1:5]
G2 <- as.data.frame(G2)
str(G2)
numo <- atcg1234(data=G2, ploidy=4); 
numo$M[1:5,1:5]
numo$ref.allele[,1:5]
have.both = intersect(rownames(numo$M), Y$all)

G3 <- numo$M[have.both,]
G3[1:5,1:5]
dim(G3)

Y1 <- Y[match(have.both,Y$all),]
dim(Y1)

G_OK <- Gmatrix(G3, method = "Slater", ploidy=4)
det(G_OK)
G_OK[1:5,1:5]
G_OK_Inv <- solve(G_OK)
PCs <- prcomp(G_OK, scale=TRUE)
head(PCs)
class(PCs)
PCA1 <- autoplot(PCs)
PCA1 + theme_bw()

autoplot(PCs)

data.3 <- merge(as.data.frame(pheno1), as.data.frame(G_OK), by = 'row.names', all = TRUE)
data.3 <- data.3 %>% remove_rownames() %>% column_to_rownames(var = 'Row.names')
data.3[1:5,1:5]
dim(data.3) # 265 6797

G_OK[1:5,1:5]

PCs <- prcomp(G_OK, scale=TRUE)
head(PCs)
class(PCs)
PCA1 <- autoplot(PCs)
PCA1 + theme_bw()
