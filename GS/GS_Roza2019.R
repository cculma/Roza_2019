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

setwd("~/Documents/git/big_files/")
Y <- read.csv("BLUE_Yi_sqrt_SpATS_DArT.csv", header = T)

colnames(Y)
head(Y)
Y <- Y[,c(1:19)]
Y1 <- na.omit(Y) 

G <- read.table("Roza2019_06_GS.txt", header = TRUE, check.names = F)
G[1:5,1:5]
dim(G) # 424 51081
common <- intersect(Y1$Roza2019_VCF_ID,rownames(G))

marks <- G[common,]
marks[1:5,1:5]
dim(marks) # 424 51081
class(marks)

marks.1 <- marks %>% dplyr::select(contains("chr1.1_"))
marks.2 <- marks %>% dplyr::select(contains("chr2.1_"))
marks.3 <- marks %>% dplyr::select(contains("chr3.1_"))
marks.4 <- marks %>% dplyr::select(contains("chr4.1_"))
marks.5 <- marks %>% dplyr::select(contains("chr5.1_"))
marks.6 <- marks %>% dplyr::select(contains("chr6.1_"))
marks.7 <- marks %>% dplyr::select(contains("chr7.1_"))
marks.8 <- marks %>% dplyr::select(contains("chr8.1_"))



Y2 <- Y1[match(common, Y1$Roza2019_VCF_ID),]
dim(Y2) # 424  19
Y2 <- Y2 %>% remove_rownames() %>% column_to_rownames(var = "Roza2019_VCF_ID")
colnames(Y2)
Y3 <- Y2[,1,drop = F]
data <- merge(as.data.frame(Y3), as.data.frame(marks.1), by = 'row.names', all = TRUE)  %>% column_to_rownames(var = 'Row.names')
data[1:5,1:5]
data <- na.omit(data)

# var_imp -----------------------------------------------------------------

ctrl <- trainControl(method = 'cv', number = 10, savePredictions = 'final')
svm1 <- train(may_20 ~., data = data, method = "svmRadial", trControl = ctrl)
imp_svm1 <- varImp(svm1, scale = T)$importance
imp_svm1 <- varImp(svm1, scale = T)
imp_svm1 <- imp_svm1$importance
hist(imp_svm1$Overall)
colnames(imp_svm1)[1] <- "svm"

colnames(pheno1)[1] <- "pheno"
print(dim(pheno1))
pheno1 <- na.omit(pheno1)

svm1 <- train(pheno ~., data = pheno1, method = "svmRadial", trControl = ctrl, tuneLength = 4)




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
