rm(list = ls())

library(GWASpoly)
library(tidyverse)
library(vcfR)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(data.table)
library(ggthemes)
library(hrbrthemes)
library(plotly)
library(GenomicRanges)
library(genomation)
library(plyranges)
library(Repitools)
library(devtools)
library(sommer)

####################
setwd("~/Documents/git/big_files/")

pheno <- read.csv("model2f.csv")
head(pheno)
trait1 <- c("BLUP")
trait1
params <- set.params(fixed=c("PC1","PC2","PC3"),
                     fixed.type= c("numeric","numeric","numeric"), n.PC = 3)

models_1 <- c("general", "additive", "1-dom", "2-dom",  "diplo-additive", "diplo-general")

data_1 <- read.GWASpoly(ploidy=4, 
                        pheno.file="model2f.csv", 
                        geno.file="Roza2019_06_GWASPoly.txt", 
                        format="numeric", n.traits=length(trait1), delim=",")


data_2 <- set.K(data = data_1, LOCO = T, n.core = 50)
data_3 <- GWASpoly(data = data_2, models = models_1, traits = trait1, params = params, n.core = 50)

# save(Yi_data_3, file = "~/Documents/git/big_files/yield_DS.RData")
# save.image("~/Documents/git/big_files/yield_DS.RData")
# load("~/Documents/git/big_files/data_Yi_DS_20.RData")

data_5 <- set.threshold(data_3, method= "Bonferroni", level=0.05)
QTL_01 <- get.QTL(data_5)

QTL_02 <- QTL_01 %>% distinct(Marker, .keep_all = T)


QTL_06.1 <- QTL_01 %>% distinct(Marker, .keep_all = T)
QTL_06.1 <- QTL_06.1 %>% unite(col = "marker1", c("Chrom", "Position"), sep = "_", remove = F)

lev0 <- unique(QTL_01$Trait)
manhattan.plot(data = data_5)



M3 <- manhattan.plot(data = data_5, traits = lev0) + theme_classic(base_family = "Arial", base_size = 12) + theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y = element_text(size = 12), plot.tag = element_blank()) 
M3

ggsave(filename = "~/Documents/git/big_files/M0_sqrt_DS.jpg", plot = M3, width = 8, height = 8)

write.csv(QTL_02, "markers_yield_DS.csv", quote = F, row.names = F)




setwd("~/Documents/git/big_files/")
Y <- pheno
colnames(Y)
head(Y)
Y <- Y %>% dplyr::select(Roza2019_VCF_ID, BLUP)
Y1 <- na.omit(Y) 
head(Y1)

G <- read.table("Roza2019_06_GS.txt", header = TRUE, check.names = F)
G[1:5,1:5]
dim(G) # 424 51081
common <- intersect(Y1$Roza2019_VCF_ID,rownames(G))

marks <- G[common,]
marks[1:5,1:5]
dim(marks) # 424 51081
class(marks)

marks.1 <- marks %>% dplyr::select(QTL_06.1$marker1)
dim(marks.1) # 424  20

Y2 <- Y1[match(common, Y1$Roza2019_VCF_ID),]
dim(Y2) # 424  2
Y2 <- Y2 %>% remove_rownames() %>% column_to_rownames(var = "Roza2019_VCF_ID")
colnames(Y2)
Y3 <- Y2

data <- merge(as.data.frame(Y3), as.data.frame(marks.1), by = 'row.names', all = TRUE) %>% column_to_rownames(var = 'Row.names')
data[1:5,1:5]
colnames(data)[1] <- "PHENOTYPE"
data <- na.omit(data)
dim(data) # 424  21

setwd("~/Documents/git/big_files/")
write.csv(data, "epistasis_model2f.csv", quote = F, row.names = T)
