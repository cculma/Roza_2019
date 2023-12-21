
library("car")
setwd("~/Documents/git/big_files/")

pheno <- read.csv("Yi_st1.csv")
pheno <- read.csv("Yi_st2.csv")
pheno <- read.csv("predictions_2.1stage_ASReml_mrbean.csv")
pheno1 <- pheno %>% dplyr::filter(trial == "aug_22") 
qqPlot(pheno1$predicted.value)

Marker <- seq(1:nrow(pheno))
pheno1 <- cbind(Marker, pheno)

?qqPlot
qqPlot(pheno$aug_22)
pheno1 <- pheno[-c(347,84,299,212,363,172,224,178,331,168),]
pheno1 <- pheno[-c(347,84,299,212,363),]
pheno1 <- pheno[-c(347,84,197,172,363,322),]
qqPlot(pheno1$aug_22)


write.csv(pheno1, "Yi_st1.1.csv", row.names = F, quote = F)


setwd("~/Documents/git/big_files/")

pheno <- read.csv("Yi_st1.csv", row.names = 1)
pheno <- read.csv("Yi_st1.1.csv", row.names = 1)
head(pheno)
trait1 <- colnames(pheno)[1:(length(colnames(pheno))-3)]
trait1
params <- set.params(fixed=c("PC1","PC2","PC3"),
                     fixed.type=rep("numeric",3), n.PC = 3)
models_1 <- c("general", "additive", "1-dom", "2-dom",  "diplo-additive", "diplo-general")

data_1 <- read.GWASpoly(ploidy=4, 
                        pheno.file="Yi_st1.1.csv", 
                        geno.file="Roza2019_06_GWASPoly.txt", 
                        format="numeric", n.traits=length(trait1), delim=",")

data_2 <- set.K(data = data_1, LOCO = T, n.core = 32)
Yi_data_3 <- GWASpoly(data = data_2, models = models_1, traits = "aug_22", params = params, n.core = 60)

data_5 <- set.threshold(Yi_data_3, method= "Bonferroni", level=0.05)
QTL_01 <- get.QTL(data_5)
manhattan.plot(data = data_5)
QTL_02 <- QTL_01 %>% distinct(Marker, .keep_all = T) 
manhattan.plot(data = data_5, traits = "aug_22", models = c("general","diplo-general","1-dom", "2-dom"))

QTL_03 <- QTL_01 %>% dplyr::filter(Trait == "aug_22") %>% dplyr::filter(Model == "general") %>% distinct(Marker, .keep_all = T) 
manhattan.plot(data = data_5, traits = NULL, models = "general")

