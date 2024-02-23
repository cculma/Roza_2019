# conver "ATGC" format to numeric with Sommer

#################
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
cl <- makePSOCKcluster(15)
registerDoParallel(cl)
#################

setwd("~/Documents/Cesar/blup_data/Roza2019/GS1")

setwd("~/Documents/git/big_files/")
# Y <- read.csv("/home/hawkins/Documents/Cesar/blup_data/Roza2019/GWASPOLY/BLUP_BLUE.1.csv", header = T)

Y <- read.csv("~/Documents/git/big_files/BLUE_Yi_sqrt_SpATS_DArT.csv", header = T)

#Y <- as.matrix(read.csv("/home/hawkins/Documents/Cesar/blup_data/Roza2019/GWASPOLY/BLUP_BLUE.1.csv", header = TRUE, row.names = 1))
colnames(Y)
head(Y)
Y <- Y[,-c(24:26)]
Y1 <- na.omit(Y) 


# G <- read.csv("/home/hawkins/Documents/Cesar/NGSEP/ngsep_tutorial/Roza_2019/2_vcf/MPP_Ms2_GWASPoly.txt", header = TRUE, check.names = F)
G <- read.table("Roza2019_06_GS.txt", header = TRUE, check.names = F)

G[1:5,1:5]
dim(G)
# G <- G[,-1]
# G1 <- unite(G, SNP, 1,2, sep = "_", remove = TRUE, na.rm = FALSE)
# G1[1:5,1:5]
# G2 <- t(G1)
# G1 <- G1 %>% remove_rownames() %>% column_to_rownames(var = 'SNP')
# str(G2)
# G2[1:5,1:5]
# rownames(G2)
# colnames(G2)
# G2 <- as.data.frame(G2)
# numo <- atcg1234(data=G2, ploidy=4)
# 
# numo$M[1:5,1:5]; 
# numo$ref.allele[,1:5]

common <- intersect(Y1$Roza2019_VCF_ID,rownames(G))

marks <- G[common,]
marks[1:5,1:5]
dim(marks)
Y2 <- Y1[match(common, Y1$Roza2019_VCF_ID),]
dim(Y2)
dim(marks)

Y2 <- Y2 %>% remove_rownames() %>% column_to_rownames(var = "Roza2019_VCF_ID")
colnames(Y2)
Y3 <- Y2[,1,drop = F]
data <- merge(as.data.frame(Y3), as.data.frame(marks), by = 'row.names', all = TRUE)
data[1:5,1:5]
Y4 <- as.matrix(data %>% remove_rownames() %>% column_to_rownames(var = 'Row.names'))
Y4[1:5,1:5]
dim(Y4)


#################

library(GROAN)
#creating a GROAN.NoisyDataset without any extra noise injected

roza1 = createNoisyDataset(
  name = 'may_20',
  genotypes = Y4[,-1], 
  phenotypes = Y4[,1],
  ploidy = 4
)
print(roza1)

wb = createWorkbench(
  folds = 10, reps = 1, stratified = FALSE, 
  outfolder = NULL, saveHyperParms = FALSE, saveExtraData = FALSE,
  regressor = phenoRegressor.rrBLUP, regressor.name = 'rrBLUP'
)
print(wb)
res = GROAN.run(roza1, wb)
sum.res <- summary(res)


GS1 <- as.matrix(marks)
G_Van <- Gmatrix(GS1, method="VanRaden", ploidy=4) 
# G_FullAutopolyploid <- Gmatrix(marks, method="Slater", ploidy=4)
G_Van[1:5,1:5]

roza2 = createNoisyDataset(
  name = 'yield',
  #  genotypes = data[,-1], 
  phenotypes = phenos,
  covariance = G_Van,
  ploidy = 4)
  
roza2 = createNoisyDataset(
  name = 'yield',
  genotypes = NULL, 
  phenotypes = Y4[,1],
  covariance = G_Van,
  ploidy = 4
)
print(roza2)
res2 = GROAN.run(roza2, wb)
sum.res2 <- summary(res2)

res = GROAN.run(roza2, wb)
print(res[,c("dataset.train", "dataset.test", "pearson")])

# actual predictions ------------------------------------------------------

#GROAN.KI has 103 samples, we'll use the first 50 samples for training
#and the remaining will be predicted
my.pheno = Y3$may_20
nrow(my.pheno)
my.pheno[200:424] = NA

#doing the predictions
res = phenoRegressor.rrBLUP(phenotypes = my.pheno, genotypes = marks)

#we just obtained a list with the following fields
print(names(res))


plot(
  x = my.pheno[200:424], xlab = "yield", 
  y = res$predictions[200:424], ylab = "Predicted values"
)
abline(a=0, b=1) 
