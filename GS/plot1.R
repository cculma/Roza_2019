rm(list = ls())
library(rrBLUP)
library(caret)
library(GROAN)
library(tidyverse)
library(hrbrthemes)
library(grid)

# load matrix 424 80137 with MET_seq in the first col
load("~/Documents/Cesar/blup_data/Roza2019/Analysis_2021/GS/data.3.RData")
rm(data.3)
# load all WGWBLUP matrix 424 424 using GWASpoly scores 8 elements for MET_sep
load("~/Documents/Cesar/blup_data/Roza2019/Analysis_2021/GS/L6_140646.RData") 

# load G3 424 140646
load("~/Documents/Cesar/blup_data/Roza2019/Analysis_2021/GS/GBLUP_matrix_140646.RData")


#######
pheno <- read.csv("~/Documents/Cesar/blup_data/Roza2019/Analysis_2021/GS/BLUP_mrbean.csv", row.names = 1)
colnames(pheno)
P1 <- pheno[,c(1,8)]
head(P1)
pheno1 <- P1[,2,drop = F]
#######

rr.model <- mixed.solve(y = pheno2, Z = NULL, K =trn, SE = FALSE, return.Hinv = TRUE)
rr.model.matrix <- as.matrix(rr.model$u)

###################
set.seed(123)
training <- createDataPartition(P1[,1], p = 0.2, list = T)
training <- createDataPartition(P1[,1], times = 3, p = 0.2, list = T)

trn <- P1[training[[1]], ]
tst <- P1[-training[[1]], ]
rownames(trn)
rownames(P1)

P2 <- cbind(a = 0, P1)
P2 <- P2[,1, drop = F]

data.3 <- merge(as.data.frame(trn), as.data.frame(P2), by = 'row.names', all = TRUE)
data.3 <- data.3 %>% remove_rownames() %>% column_to_rownames(var = 'Row.names')
data.4 <- merge(as.data.frame(tst), as.data.frame(data.3), by = 'row.names', all = TRUE)

tst <- P1[-training,]
tst1 <- L6[[2]][-training,]
is.na[P1]

rownames(tst)
rownames(tst1)

dim(trn)
dim(tst)
class(trn)


###############
# plot best models using GROAN
my.pheno = P1$MET_sep
my.pheno[324:424] = NA

res.1 = phenoRegressor.rrBLUP(phenotypes = my.pheno, genotypes = NULL, covariances = L6[[2]])
res.2 = phenoRegressor.rrBLUP(phenotypes = my.pheno, genotypes = NULL, covariances = G_3.1)
res.3 = phenoRegressor.rrBLUP(phenotypes = my.pheno, genotypes = NULL, covariances = L6[[1]])

# generate dataframe to plot scatter plot
P3 <- data.frame(P1$MET_sep[324:424], res.1$predictions[324:424])
P3 <- data.frame(P1$MET_sep[324:424], res.2$predictions[324:424])
P3 <- data.frame(P1$MET_sep[324:424], res.3$predictions[324:424])

colnames(P3) <- c('EVB', 'GEVB')
cor(P3$EVB, P3$GEVB)
P3$pc <- predict(prcomp(~EVB+GEVB, P3))[,1]


ggplot(P3, aes(x = EVB, y = GEVB, color = pc)) + geom_point(size = 3, alpha = 0.5) + theme_minimal(base_size = 16) + scale_color_gradient(low = "yellowgreen", high = "royalblue4") + theme_ipsum() + labs(title = 'Model = WGBLUP additive\n140,646 predictors used for regression', x = 'Original phenotypes\nMET_Sep BLUPs', y = 'Predicted phenotypes\nGEBV')  + annotation_custom(grob1) + theme(legend.position = 'none')

grob1 = grobTree(textGrob(paste("Pearson Correlation : ", round(cor(P3$EVB, P3$GEVB), 4) ), x = 0, y = 0.97, hjust = 0, gp = gpar(fontsize = 11)))


#############
# WGBLUP Additive using all pheno values
res.0 = phenoRegressor.rrBLUP(phenotypes = P1$MET_sep, genotypes = NULL, covariances = L6[[2]])

P0 <- data.frame(P1$MET_sep, res.0$predictions)
colnames(P0) <- c('EVB', 'GEVB')
cor(P0$EVB, P0$GEVB)

plot(GROAN.KI$yield, pch=20, main = 'True (black) vs. Noisy (red)', xlab = 'Samples', ylab = 'Phenotypes')
#plotting an instance of the phenotypes with noise injected 
points(getNoisyPhenotype(nds.normal_noise), col='red')

# what happen if I compare with other phenotypic data

res.4 = phenoRegressor.rrBLUP(phenotypes = P1$MET_all, genotypes = NULL, covariances = L6[[2]])
P0 <- data.frame(P1$MET_all, res.0$predictions)
colnames(P0) <- c('EVB', 'GEVB')
cor(P0$EVB, P0$GEVB)
P0$pc <- predict(prcomp(~EVB+GEVB, P0))[,1]

ggplot(P0, aes(x = EVB, y = GEVB, color = pc)) + geom_point(size = 3, alpha = 0.5) + theme_minimal(base_size = 16) + scale_color_gradient(low = "yellowgreen", high = "royalblue4") + theme_ipsum() + labs(title = 'Model = WGBLUP additive\n140,646 predictors used for regression', x = 'Original phenotypes\nMET_all BLUPs', y = 'Predicted phenotypes\nGEBV')  + annotation_custom(grob1) + theme(legend.position = 'none')

grob1 = grobTree(textGrob(paste("Pearson Correlation : ", round(cor(P0$EVB, P0$GEVB), 4) ), x = 0, y = 0.97, hjust = 0, gp = gpar(fontsize = 11)))

# comparison between sep and all
ggplot(P1, aes(x = MET_all, y = MET_sep, color = pc)) + geom_point(size = 3, alpha = 0.5) + theme_minimal(base_size = 16) + scale_color_gradient(low = "yellowgreen", high = "royalblue4") + theme_ipsum() + labs(title = 'Correlation of MET BLUP values', x = 'Original phenotypes\nMET_all BLUPs', y = 'Original phenotypes\nMET_sep BLUPs')  + annotation_custom(grob1) + theme(legend.position = 'none')

P1$pc <- predict(prcomp(~MET_all+MET_sep, P1))[,1]
grob1 = grobTree(textGrob(paste("Pearson Correlation : ", round(cor(P1$MET_sep, P1$MET_all), 4) ), x = 0, y = 0.97, hjust = 0, gp = gpar(fontsize = 11)))

# correlation among traits
colnames(pheno)
P2 <- pheno[1:(length(pheno)-3)]
P2.1 <- cor(P2)
P2 <- pheno[,c(1:8)]

trait1 <- colnames(pheno)[1:(length(colnames(pheno))-3)]

#################
data.3 <- merge(as.data.frame(pheno1), as.data.frame(G3), by = 'row.names', all = TRUE)
data.3 <- data.3 %>% remove_rownames() %>% column_to_rownames(var = 'Row.names')
dim(data.3) # 424 140647
data.3[1:5,1:5]
# RRBLUP
MPP_1 = createNoisyDataset(name = 'MET_sep', 
                           genotypes = data.3[,-1], 
                           phenotypes = data.3[,1], ploidy = 4)
wb_1 = createWorkbench(folds = 10, reps = 10, stratified = FALSE, outfolder = NULL, saveHyperParms = FALSE, saveExtraData = FALSE, regressor = phenoRegressor.rrBLUP, regressor.name = "RRBLUP_1")
res_1 = GROAN.run(MPP_1, wb_1)
write.csv(res_1, "~/Documents/Cesar/blup_data/Roza2019/Analysis_2021/GS/RRBLUP_140647_MET_sep.csv")

