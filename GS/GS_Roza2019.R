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
cl <- makePSOCKcluster(40)
registerDoParallel(cl)
#################

setwd("~/Documents/Cesar/blup_data/Roza2019/GS1")

setwd("~/Documents/git/big_files/")
# Y <- read.csv("/home/hawkins/Documents/Cesar/blup_data/Roza2019/GWASPOLY/BLUP_BLUE.1.csv", header = T)

Y <- read.csv("BLUE_Yi_sqrt_SpATS_DArT.csv", header = T)

colnames(Y)
head(Y)
Y <- Y[,-c(24:26)]
Y1 <- na.omit(Y) 

G <- read.table("Roza2019_06_GS.txt", header = TRUE, check.names = F)

G[1:5,1:5]
dim(G)

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

roza1 = createNoisyDataset(name = 'may_20',
  genotypes = Y4[,-1], 
  phenotypes = Y4[,1],
  ploidy = 4
)
print(roza1)

set.seed(123)
wb = createWorkbench(folds = 10, reps = 10, stratified = FALSE, outfolder = NULL, saveHyperParms = T, saveExtraData = T, regressor = phenoRegressor.rrBLUP, regressor.name = "GBLUP")

?createWorkbench
# wb = createWorkbench(folds = 10, reps = 1, stratified = FALSE, outfolder = NULL, saveHyperParms = FALSE, saveExtraData = FALSE,regressor = phenoRegressor.rrBLUP, regressor.name = 'rrBLUP')

print(wb)
res = GROAN.run(roza1, wb)
sum.res <- summary(res)


GS1 <- as.matrix(marks)
G_Van <- Gmatrix(GS1, method="VanRaden", ploidy=4) 
G_Van[1:5,1:5]


rownames(G_Van) == colnames(G_Van)
rownames(G_Van) == rownames(Y2)

head(Y2)

roza2 = createNoisyDataset(
  name = 'yield',
  genotypes = NULL, 
  phenotypes = Y2[,1],
  covariance = G_Van,
  ploidy = 4
)
print(roza2)
res2 = GROAN.run(roza2, wb)

mean(res2$pearson)
# 0.1964714

print(res[,c("dataset.train", "dataset.test", "pearson")])

# actual predictions ------------------------------------------------------

res.4 = phenoRegressor.rrBLUP(phenotypes =  Y2[,1], genotypes = NULL, covariances = G_Van, SE = T, return.Hinv = T)

P0 <- data.frame(Y2[,1], res.4$predictions)
colnames(P0) <- c('yield', 'GEVB')
P0$pc <- predict(prcomp(~yield+GEVB, P0))[,1]
cor(P0$yield, P0$GEVB)

P0 <- as.data.frame(res.4$predictions) %>% rownames_to_column("INDIV")
colnames(P0)[2] <- "GROAN"

P1 <- as.data.frame(res.4$predictions) %>% rownames_to_column("INDIV")
colnames(P1)[2] <- "GROAN1"

GEBV_ASReml <- read.csv("GEBV_GBLUP_st1.csv")
GEBV_ASReml <- GEBV_ASReml %>% dplyr::select(INDIV, may_20)
colnames(GEBV_ASReml)[2] <- "ASReml"

head(Y2)
Y3 <- Y2 %>% rownames_to_column("INDIV") %>% dplyr::select(INDIV, may_20)

Y4 <- inner_join(P0, GEBV_ASReml, by = "INDIV") %>% inner_join(., Y3, by = "INDIV") %>% inner_join(., P1, by = "INDIV") %>% column_to_rownames("INDIV") 

cor(Y4)
Y4 <- inner_join(P0, P1, by = "INDIV") 


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
