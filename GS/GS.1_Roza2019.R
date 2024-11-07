rm(list=ls()) 

library(neuralnet)
library(superpc)
library(xgboost)
library(kernlab)
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
library(RColorBrewer)
library(GROAN)
library(data.table)
library(GWASpoly)
library(AGHmatrix)
cl <- makePSOCKcluster(15)
registerDoParallel(cl)

setwd("/home/hawkins/Documents/Cesar/blup_data/Roza2019/GS1")
load("/home/hawkins/Documents/Cesar/blup_data/Roza2019/GS1/Roza2019_num_mat.RData")
load("/home/hawkins/Documents/Cesar/blup_data/Roza2019/GWASPOLY/single_hw_roza2019.RData")

dim(Roza2019_num_mat)
t26 <- Roza2019_num_mat
t26[1:5,1:10]
?write.GWASpoly
write.GWASpoly(data_6, )

LD.plot(data_6, max.pair = 10000, dof = 8)

############
# to write GWAS table scores
data_6 <- set.threshold(single_hw_roza2019, method= "Bonferroni", level=0.05)
traits <- c("MET_yield")
for (i in traits) {
  write.GWASpoly(data_6, filename=paste("MPP",i, sep="_"), 
                 trait=i, what="scores", delim = "\t")
}
MET_yield <- read.table("/home/hawkins/Documents/Cesar/blup_data/Roza2019/GS1/MPP_MET_yield", header = T, check.names = F)
# write.GWASpoly(data_6, trait = 'MET_yield', filename = 'MET_scores.csv', what = "effects", delim = ",")

#############

t4 <- c('MET_yield')
t5 <- c('0.5', '1.0', '1.5', '2.0', '2.5', '3.0', '3.5', '4.0')

# this for loop allows to filter MET_yield according to different columns with different tresholds 
t1 <- list(MET_yield)
t2 <- list()
for (i in t1) {
  for (j in seq(0.5,4,0.5)) {
    k <- i %>% filter(`general`>=j | `additive`>=j | 
                      `diplo-additive`>=j | `diplo-general`>=j | `1-dom-alt`>=j |
                      `1-dom-ref`>=j | `2-dom-alt`>=j | `2-dom-ref`>=j)
    k.1 <- k %>% unite(Chrom1, 2:3, remove = T)
    k.1 <- k.1[,2,drop = F]
    t2[[length(t2)+1]] = k.1
  }
}

my_vector = c() # to create a vector with a combination of t4 and t5. To rename the lists
for (i in t4) {
  for (j in t5) {
    t6 <- paste(i,j, sep = "_", collapse = NULL)
    my_vector=c(my_vector,t6)
  }
}
names(t2) <- my_vector
MET_yield_0.0 <- as.data.frame(colnames(Roza_G))
colnames(MET_yield_0.0) <- "Chrom1"

t3 <- list()
for (i in 1:length(t2)) {
  A1 <- t26[ ,which((colnames(t26) %in% t2[[i]] $Chrom)==TRUE)]
  t3[[length(t3)+1]] = A1
}
names(t3) <- names(t2)
t4 <- c(MET_yield_0.0 = list(t26), t3)
# t4 contains all geno matrices 
names(t4)
# load pheno data
pheno <- read.csv("/home/hawkins/Documents/Cesar/blup_data/Roza2019/GWASPOLY/Roza2019_HW.1.csv", row.names = 1)
head(pheno)
P1 <- pheno[,1, drop = F]
row.names(t4[[5]])

# join pheno data to geno data
t5 <- list()
for (i in 1:length(t4)) {
  data <- merge(P1, t4[[i]], by = 'row.names', all = FALSE)
  data <- data %>% remove_rownames() %>% column_to_rownames(var = 'Row.names')
  data <- as.matrix(data)
  t5[[length(t5)+1]] = data
}
names(t5) <- names(t4)

# save list of MET_yield with dif matrices 
MET_GS <- t5
rm(MET_GS)
save(MET_GS, file="/home/hawkins/Documents/Cesar/blup_data/Roza2019/GS1/MET_GS.RData")
load("/home/hawkins/Documents/Cesar/blup_data/Roza2019/GS1/MET_GS.RData")
t5 <- MET_GS
t5[[1]][1:5,1:5]
#################
# rrBLUP
t6 <- list()
for (i in 1:length(t5)) {
  nds = createNoisyDataset(
    name = 'MET_yield_all',
    genotypes = t5[[i]][,-1], 
    phenotypes = t5[[i]][,1],
    ploidy = 4)
  set.seed(123)
  res = GROAN.run(nds, wb)
  t6[[length(t6)+1]] = res
}
names(t6) <- names(t5)
t6 <- rbindlist(t6, use.names = T, fill = T, idcol = 'threshold')

wb = createWorkbench(
  folds = 10, reps = 1, stratified = FALSE, 
  outfolder = NULL, saveHyperParms = FALSE, saveExtraData = FALSE,
  regressor = phenoRegressor.rrBLUP, regressor.name = 'rrBLUP')
#############
# gBLUP using three matrices
# G_Van <-Gmatrix(B1,method="VanRaden",ploidy=4)
# G_Dom <-Gmatrix(B1,method="Endelman",ploidy=4)
# G_Ful <-Gmatrix(B1,method="Slater",ploidy=4)

# comment and uncomment each line acording to the matrix
t8 <- list()
for (i in 1:length(t5)) {
  data1 <- t5[[i]][,-1]
#  G_mat <-Gmatrix(data1,method="VanRaden",ploidy=4)
#  G_mat <-Gmatrix(data1,method="Endelman",ploidy=4)
  G_mat <-Gmatrix(data1,method="Slater",ploidy=4)
  t8[[length(t8)+1]] = G_mat
}
names(t8) <- names(t5)

t8.1 <- list()
for (i in 1:length(t8)) {
  data1 <- t8[[i]][,-1]
  G_mat <-Gmatrix(data1,method="Slater",ploidy=4)
  t8.1[[length(t8.1)+1]] = G_mat
}
names(t8.1) <- names(t5)
# matrix for phenotypic trait 'MET_yield'
y <- t5[[1]][,1]
y <- matrix(y, ncol = 1)
colnames(y) <- 'MET_yield'
rownames(y) <- rownames(t5[[1]])

t9 <- list()
for (i in 1:length(t8.1)) {
  G.BLUP = createNoisyDataset(
    name = 'MET_yield_additive',
    phenotypes = y,
    genotypes = NULL,
#    covariance = t8[[i]],
    covariance = t8.1[[i]],
    ploidy = 4)
  set.seed(123)
  res.2 = GROAN.run(G.BLUP, wb)
  t9[[length(t9)+1]] = res.2
}
names(t9) <- names(t8)
t9 <- rbindlist(t9, use.names = T, fill = T, idcol = 'threshold')
# t7 <- rbind(t9, t7)

write.csv(t7, file = "/home/hawkins/Documents/Cesar/blup_data/Roza2019/GS1/GBLUP_MET.csv", quote = F, row.names = F)

##################
# change Workbench according to G_matrix model
# wb = createWorkbench(
#   folds = 10, reps = 1, stratified = FALSE, 
#   outfolder = NULL, saveHyperParms = FALSE, saveExtraData = FALSE,
#   regressor = phenoRegressor.rrBLUP, regressor.name = 'gBLUP_VanRaden'
# ) # VanRaden
# wb = createWorkbench(
#   folds = 10, reps = 1, stratified = FALSE, 
#   outfolder = NULL, saveHyperParms = FALSE, saveExtraData = FALSE,
#   regressor = phenoRegressor.rrBLUP, regressor.name = 'gBLUP_Dominance'
# ) # Endelman
wb = createWorkbench(
  folds = 10, reps = 1, stratified = FALSE, 
  outfolder = NULL, saveHyperParms = FALSE, saveExtraData = FALSE,
  regressor = phenoRegressor.rrBLUP, regressor.name = 'gBLUP_FullAuto'
) # "Slater"

# plot hist y
y <- as.data.frame(y)
ggplot(y, aes(x=MET_yield)) + geom_histogram(color="black", fill="yellowgreen") + theme_classic(base_size = 16)


#############
# accuracy plot 
t7 <- read.csv( "/home/hawkins/Documents/Cesar/blup_data/Roza2019/GS1/GBLUP_MET.csv")
t8 <- separate(t7, 1, c('trait_1', 'trait_2', 'p_value'), sep = '_', remove = TRUE, convert = FALSE, extra = 'warn')
cnames <- c('p_value', 'regressor')
t8[cnames] <- lapply(t8[cnames], factor)

t8$p_value <- as.factor(t8$p_value)
t8$regressor <- as.factor(t8$regressor)
str(t8)

ggplot(data = t8, aes(x=p_value, y = pearson, group = regressor)) + geom_line(aes(color = regressor)) + geom_point(aes(color = regressor)) + theme_classic(base_size = 16) + scale_color_brewer(palette = "Dark2") 

#############
# to correlate p_values according to differnet model used (GWASpoly) 
str(MET_yield)
k <- MET_yield %>% unite(Chrom1, 2:3, remove = T)
colnames(k)
k <- k[,c(2,5:12)]
k <- k %>% remove_rownames() %>% column_to_rownames(var = 'Chrom1')
k <- as.matrix(k)
class(k)
str(k)
cor(k)
library(psych)
pairs.panels(k)
library(PerformanceAnalytics)
chart.Correlation(k)

##############
# to evaluate if there are differences using different models

# k <- i %>% filter(`general`>=j | `additive`>=j |
#                   `diplo-additive`>=j | `diplo-general`>=j | `1-dom-alt`>=j |
#                   `1-dom-ref`>=j | `2-dom-alt`>=j | `2-dom-ref`>=j)
colnames(MET_yield)
t1 <- list(MET_yield)
t8 <- list()
for (i in t1) {
  for (j in seq(0.5,4,0.5)) {
    k <- i %>% filter(`additive`>=j)
    k.1 <- k %>% unite(Chrom1, 2:3, remove = T)
    k.1 <- k.1[,2,drop = F]
    A1 <- t26[ ,which((colnames(t26) %in% k.1$Chrom1)==TRUE)]
    data <- merge(P1, A1, by = 'row.names', all = FALSE)
    data <- data %>% remove_rownames() %>% column_to_rownames(var = 'Row.names')
    data <- as.matrix(data)
    t8[[length(t8)+1]] = data
  }
}
names(t8) <- names(t2)
data <- merge(P1, t26, by = 'row.names', all = FALSE)
data <- data %>% remove_rownames() %>% column_to_rownames(var = 'Row.names')
data <- as.matrix(data)
t8 <- c(MET_yield_0.0 = list(data), t8)

# RRBLUP different models
t4 <- list()
for (i in 1:length(t8)) {
  nds = createNoisyDataset(
    name = 'MET_yield_additive',
    genotypes = t8[[i]][,-1], 
    phenotypes = t8[[i]][,1],
    ploidy = 4)
  set.seed(123)
  res = GROAN.run(nds, wb)
  t4[[length(t4)+1]] = res
}
names(t4) <- names(t8)
t4 <- rbindlist(t4, use.names = T, fill = T, idcol = 'threshold')

# t3 <- rbind(t4, t6)
t3 <- rbind(t9, t3)

write.csv(t3, file = "/home/hawkins/Documents/Cesar/blup_data/Roza2019/GS1/rrBLUP_MET_models.csv", quote = F, row.names = F)

write.csv(t3, file = "/home/hawkins/Documents/Cesar/blup_data/Roza2019/GS1/GBLUP_MET_models.csv", quote = F, row.names = F)

wb = createWorkbench(
  folds = 10, reps = 1, stratified = FALSE, 
  outfolder = NULL, saveHyperParms = FALSE, saveExtraData = FALSE,
  regressor = phenoRegressor.rrBLUP, regressor.name = 'rrBLUP')

colnames(t3)
t4 <- separate(t3, 1, c('trait_1', 'trait_2', 'p_value'), sep = '_', remove = TRUE, convert = FALSE, extra = 'warn')
str(t4)
head(t4)
t4$dataset.train <- gsub('MET_yield.', '', t4$dataset.train)
colnames(t4)
t4 <- t4 %>% unite(model, 4:5, remove = T)

# cnames <- c('p_value', 'regressor')
# t4[cnames] <- lapply(t4[cnames], factor)
t4$p_value <- as.factor(t4$p_value)
t4$model <- as.factor(t4$model)
t4$dataset.train <- as.factor(t4$dataset.train)
head(t4)

colnames(t4)
colnames(t4)[5] <- "model"

ggplot(data = t4, aes(x=p_value, y = pearson, group = model)) + geom_line(aes(color = model)) + geom_point(aes(color = model)) + theme_classic(base_size = 16) +  scale_color_brewer(palette = "Dark2") 

library(RColorBrewer)
nb.colors <- 9
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.colors)
scale_fill_manual(values = mycolors)
#############

M <- as.matrix(A1)
G_Ful <-Gmatrix(A1,method="Slater",ploidy=4)

control <- rfeControl(functions = rfFuncs, method = 'cv', number = 10)
results <- rfe(A1[,-1], A1[,1], sizes = c(10, 500, 2000, 4000, 6000, 8000, 10000), rfeControl = control)
print(results)
predictors(results)
plot(results, type=c('g','o'))

results[["pred"]]$pred

pearson = cor(model.ranger[['pred']]$pred, model.ranger[['pred']]$obs)

ctrl <- trainControl(method = "cv", savePred=T, number=3)
model.ranger <- train(BLUE ~., data = A1, method = 'lars', trControl = train.control, use.Gram=FALSE)

cvresult <- train(BLUE~.,
                  data=A1,
                  method = "lars",
                  trControl = ctrl,
                  metric="RMSE")

coeffs <- predict.lars(cvresult$finalModel,type="coefficients")
models <- as.data.table(coeffs$coefficients)
winnermodelscoeffs <- models[which(coeffs$fraction==cvresult$bestTune$fraction)]

apply(models, 2, var)
which(apply(models, 2, var) == 0)
models.1 <- models[ - as.numeric(which(apply(models, 2, var) == 0))]

################
dim(A1)
geno <- A1[,-1]
pheno <- A1[,1,drop = F]
#training-60% validation-40%
train1 <- as.matrix(sample(1:199, 80)) #200*0.4
test1 <- setdiff(1:199,train1)
Pheno_train1 <- pheno[train1,]
m_train1 <- geno[train1,]
Pheno_valid1 <- pheno[test1,]
m_valid1 <- geno[test1,]
dim(m_valid1)

WA_Jul_18_answer <- mixed.solve(Pheno_train1, Z = m_train1, K = NULL, SE = FALSE,
                                return.Hinv = FALSE)

WA_Jul_18 <- (Pheno_train1[,1])
WA_Jul_18_answer <- mixed.solve(WA_Jul_18, Z = m_train1, K = NULL, SE = FALSE,
                                return.Hinv = FALSE)

WA_Jul_18_U <- WA_Jul_18_answer$u
WA_Jul_18_matrix <- as.matrix(WA_Jul_18_U)
dim(WA_Jul_18_matrix)
pred_WA_Jul_18_valid <- (m_valid1 %*% WA_Jul_18_matrix)
pred_WA_Jul_18 <- (pred_WA_Jul_18_valid[,1])+WA_Jul_18_answer$beta
pred_WA_Jul_18
WA_Jul_18_valid <- Pheno_valid1[,1]
WA_Jul_18_accuracy <- cor(pred_WA_Jul_18_valid, WA_Jul_18_valid, use = "complete")
WA_Jul_18_accuracy
#[1,] 0.08183855

################

traits=1
cycles=500
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  train1 <- as.matrix(sample(1:200, 80)) #200*0.4
  test1 <- setdiff(1:200,train1)
  Pheno_train1 <- Pheno_WA[train1,]
  m_train1 <- Markers_WA[train1,]
  Pheno_valid1 <- Pheno_WA[test1,]
  m_valid1 <- Markers_WA[test1,]
  
  WA_Jul_18 <- (Pheno_train1[,1])
  WA_Jul_18_answer <- mixed.solve(WA_Jul_18, Z = m_train1, K = NULL, SE = FALSE,
                                  return.Hinv = FALSE)
  WA_Jul_18_U <- WA_Jul_18_answer$u
  WA_Jul_18_matrix <- as.matrix(WA_Jul_18_U)
  pred_WA_Jul_18_valid <- (m_valid1 %*% WA_Jul_18_matrix)
  #pred_WA_Jul_18 <- (pred_WA_Jul_18_valid[,1])+WA_Jul_18_answer$beta
  #pred_WA_Jul_18
  WA_Jul_18_valid <- Pheno_valid1[,1]
  WA_Jul_18_accuracy[r,1] <- cor(pred_WA_Jul_18_valid, WA_Jul_18_valid, use = "complete")
  
}
mean(WA_Jul_18_accuracy)
install.packages("synbreed")
install.packages('/home/hawkins/Downloads/synbreed_0.12-14.tar.gz',repos = NULL)

library(Matrix)
library(synbreed)

data <- rnorm(1e6)
zero_index <- sample(1e6)[1:9e5]
data[zero_index] <- 0
mat <- matrix(data, ncol = 1000)
mat[1:5,1:5]
print(object.size(mat),units="auto")
mat_sparse <- Matrix(mat, sparse=TRUE)
mat_sparse[1:5,1:5]

dim(geno)
is.matrix(geno)
geno.1 <- as.matrix(geno)
geno.1[1:5,1:5]


#######

Y1 <- read.csv("/home/hawkins/Documents/Cesar/blup_data/Roza2019/yield/Roza_2019_yield_2020.csv")
Y2 <- read.csv("/home/hawkins/Documents/Cesar/blup_data/Roza2019/hw/Roza_2019_HW.csv")
colnames(Y1)
colnames(Y2)
Y1 <- Y1[,c(2,4,5,6)]
Y3 <- inner_join(Y1, Y2, by = "acc")
colnames(Y3)
Y3 <- Y3[,-4]

write.csv(Y3, "/home/hawkins/Documents/Cesar/blup_data/Roza2019/hw/Roza_2019_HW.1.csv", row.names = F, quote = F)

str(Y3)
cols <- c('col', 'row', 'acc_num', 'rep', 'plant_num')
Y3[cols] <- lapply(Y3[cols], factor)

Y3.1 <- subset(Y3, rep == "1")
Y3.2 <- subset(Y3, rep == "2")
Y3.3 <- subset(Y3, rep == "3")

cols <- c('acc', 'col', 'row','rep')
Y1[cols] <- lapply(Y1[cols], factor)
str(Y1)
Y1.1 <- subset(Y1, rep == "1")
str(Y1.1)
library(data.table) 
A <- setDT(Y3.1)[factor == "acc_num", mean(width_2020)]
## [1] 1.5
str(Y3.1)
A <- tapply(Y3.1$height_2019, Y3.1$acc, mean)
length(A)
tapply(df$speed, df$dive, mean)

A <- read.csv("/home/hawkins/Documents/Cesar/blup_data/Roza2019/yield/MET_Roza_2019_yield.csv")
colnames(A)
str(A)
A$gen <- as.factor(A$gen)
B <- read.csv("/home/hawkins/Documents/Cesar/blup_data/Roza2019/GWASPOLY/id_gbs.csv")
colnames(B) <- c("all", "gen")
str(B)
C <- inner_join(B, A, by = "gen")
colnames(C)
C <- C[,c(1,3)]
colnames(C) <- c("all", "MET_yield")
D <- read.csv("/home/hawkins/Documents/Cesar/blup_data/Roza2019/GWASPOLY/Roza2019_HW.csv")
E <- inner_join(C, D, by = "all")

write.csv(E, "/home/hawkins/Documents/Cesar/blup_data/Roza2019/GWASPOLY/Roza2019_HW.1.csv", row.names = F, quote = F)


###
# start GS analysis

G <- read.csv("/home/hawkins/Documents/Cesar/blup_data/Roza2019/GWASPOLY/MPP_Ms2_GWASPoly.txt", header = TRUE, row.names = 1, check.names = F)

G[1:5,1:5]
G1 <- G %>% unite(Chrom1, 1:2, remove = T)
G1[1:5,1:5]
G3 <- as.matrix(G1 %>% remove_rownames() %>% column_to_rownames(var = "Chrom1"))
G3 <- t(G3)
G3[1:5,1:5]
G3 <- as.data.frame(G3)
numo.1 <- atcg1234(data=G3, ploidy=4)
G4 <- numo.1$M
G4[1:5,1:5]

dim(G4)
Roza2019_num_mat <- G4
save(Roza2019_num_mat, file = "/home/hawkins/Documents/Cesar/blup_data/Roza2019/GS1/Roza2019_num_mat.RData")


G5 <- numo.1$ref.alleles
numo.1$ref.alleles[1:5]

P1 <- read.csv("/home/hawkins/Documents/Cesar/blup_data/Roza2019/GWASPOLY/Roza2019_HW.1.csv")

head(P1)
P1 <- na.omit(P1)
P2 <- as.matrix(P1 %>% remove_rownames() %>% column_to_rownames(var = "all"))
head(P2)
str(P2)
data <- merge(as.data.frame(P2[,7]), as.data.frame(G4), by = 'row.names', all = FALSE)
data <- data %>% remove_rownames() %>% column_to_rownames(var = 'Row.names')
data[1:5,1:5]
#str(data)

# to run gblup according to book chapter 11
# this is to calculate the gmatrix from scratch

P3 <- read.table("/home/hawkins/Documents/Cesar/blup_data/Roza2019/GWASPOLY/MPP_Ms2_stats.1.txt", sep = '\t', header = F)
head(P3)
colnames(P3) <- c('chrom', 'snp1', 'snp2', 'ref_allele', 'output')
P4 <- P3[-c(3,4)]
head(P4)
P5 <- separate(P4, 3, c('n_sample', 'EH', 'OH', '1_OH_EH', 'MAF', 'RAF', 'X2_HWE', 'p_val_X2'), sep = ":", remove = TRUE, convert = FALSE, extra = "warn")
head(P5)
str(P5)
P5$MAF <- as.numeric(P5$MAF)
P6 <- P5[,c(1,2,7)]
head(P6)
dim(P6)

######
# here I use 2 as formula.
maf <- P5$MAF
str(maf)
P <- 2*rep(1, 499)%*%t(maf)
dim(P)
denom = as.numeric(2*t(maf)%*%(1 - maf))
Z <- G4-P
G = (Z%*%t(Z))/denom
dim(G)
G[1:5,1:5]

######
# here I change to 4 
P <- 4*rep(1, 499)%*%t(maf)
dim(P)
denom = as.numeric(4*t(maf)%*%(1 - maf))
Z <- G4-P
G1 = (Z%*%t(Z))/denom
dim(G1)
G1[1:5,1:5]



maf <- c(0.01, 0.1, 0.25, 0.48)
P <- 2*rep(1,3)%*%t(maf)
P2 <- 4*rep(1,3)%*%t(maf)
M <- matrix(c(0,1,0,2,2,1,1,1,2,0,0,0), nrow=3,ncol=4, byrow=T)
Z <- M-P
Z2 <- M-P2

denom = as.numeric(2*t(maf)%*%(1 - maf))
G = (Z%*%t(Z))/denom
G2 = (Z2%*%t(Z2))/denom



library(synbreed)
library(asreml)
library(GeneticsPed)

Gtab = write.relationshipMatrix(G, sorting = "ASReml", type = "none")
colnames(Gtab)[3] = "G"
head(Gtab)
G[1:5,1:5]
str(Gtab)

G.inv <- write.relationshipMatrix(G, file = NULL, sorting = c('ASReml'), type = c('ginv'), digits = 10)
head(attr(G.inv, 'rowNames'))
names(G.inv) <- c('row', 'column', 'coefficient')
head(G.inv)



dereg<-function(gs,ps,gm,pm,gi,pi,lambda,c){
  ###############################################################
  # calculates deregressed BVs and weights starting from EBV and reliabilities
  # gs, gm, gi are the EBVs of sire dam and individual respectively
  # ps, pm and pi are the reliabilities of sire, dam and individual respectively
  #lambda is the Ve/Va ratio (obtained from the BLUP analysis)
  # c is the proportion of variance (un)explained by the markers
  # returns a list of de-regressed BV accuracies of de-regressed and weights
  ################################################################
  rpa<-(ps+pm)/4
  gpa<-(gs+gm)/2
  alpha<-1/(0.5-rpa)
  delta<-(0.5-rpa)/(1-pi)
  ZZpa<-lambda*(0.5*alpha-4)+0.5*lambda*sqrt(alpha^2+16/delta)
  ZZi<-delta*ZZpa+2*lambda*(2*delta-1)
  LHS<-rbind(cbind(ZZpa+4*lambda,-2*lambda),cbind(-
                                                    2*lambda,ZZi+2*lambda))
  RHS<-rbind(gpa,gi)
  so<-LHS%*%RHS
  de<-so[2]/ZZi
  rw<-1-lambda/(ZZi+lambda)
  w<-(1/(c+(1-rw)/rw))*lambda
  return(c(de,rw,w))
}

Dominance <- kin(gp.num, ret="dom")
summary(Dominance)


class(gp.num)
gp.num$pheno[1:4,1:4]

load("/home/hawkins/Documents/Cesar/blup_data/Roza2019/GS1/Roza2019_num_mat.RData")
Roza_G <- Roza2019_num_mat
Roza_G[1:5,1:5]
is.matrix(Roza_G)
dim(Roza_G)

Roza_P <- read.csv("/home/hawkins/Documents/Cesar/blup_data/Roza2019/GWASPOLY/Roza2019_HW.1.csv", row.names = 1)
is.data.frame(Roza_P)
colnames(Roza_P)
Roza_P1 <- Roza_P[,c(1:7)]
colnames(Roza_P1)

# create.gp.Data to create genomic prediction data object
gp.roza <- create.gpData(pheno = Roza_P, geno = Roza_G, map = NULL, pedigree = NULL,
                         family = NULL, covar = NULL, reorderMap = TRUE, map.unit = "cM",
                         repeated  = NULL, modCovar = NULL, na.string="NA", cores=1)
class(gp.roza)

data <- merge(as.data.frame(Roza_P1), as.data.frame(Roza_G), by = 'row.names', all = FALSE)
data <- data %>% remove_rownames() %>% column_to_rownames(var = 'Row.names')
data[1:5,1:10]
Roza_G1 <- as.matrix(data[,-c(1:7)])
Roza_P1 <- data[,c(1:7)]
class(Roza_G1)
class(Roza_P1)
dim(Roza_G1)


G_Van <-Gmatrix(Roza_G1,method="VanRaden",ploidy=4)
G_Dom <-Gmatrix(Roza_G1,method="Endelman",ploidy=4)
G_Ful <-Gmatrix(Roza_G1,method="Slater",ploidy=4)

dim(G_Van)

summary(cv.Ablup)
summary(cv.Gblup)
class(cv.Ablup)

colnames(Roza_P1)
G <- G_Van # when Ginv <- solve(G) Error in solve.default(G) : system is computationally singular: reciprocal condition number = 1.9618e-16

G <- G_Dom
# G <- G_Ful
Ginv <- solve(G)
Ginv[1:5,1:5]
Gdia <- G%*%Ginv
Gdia[1:5,1:5]

?write.relationshipMatrix
G1 = write.relationshipMatrix(G, sorting = "ASReml", type = "ginv", digits = 2)
head(G1)

y <- Roza_P1$MET_yield
rownames(G)
rownames(Roza_P1)
dim(G)

dim(Roza_P)
head(Roza_P)

modelGBLUP<-asreml(fixed=JUI_MOT~1,
                   random=~vm(INDIV,ginv),
                   workspace=128e06,
                   data=datag)
help(vm)
