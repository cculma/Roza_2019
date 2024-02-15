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


pheno <- read.csv("BLUE_Yi_sqrt_DArT_DS1.csv")
pheno$X20
head(pheno)
trait1 <- c("X20","X21","X22","X23")
trait1 <- colnames(pheno)[2:(length(colnames(pheno))-5)]
trait1

params <- set.params(fixed=c("PC1","PC2","PC3","Stress"),
                     fixed.type= c("numeric","numeric","numeric","factor"), n.PC = 3)

models_1 <- c("general", "additive", "1-dom", "2-dom",  "diplo-additive", "diplo-general")


data_1 <- read.GWASpoly(ploidy=4, 
                        pheno.file="BLUE_Yi_sqrt_DArT_DS1.csv", 
                        geno.file="Roza2019_06_GWASPoly.txt", 
                        format="numeric", n.traits=length(trait1), delim=",")

data_2 <- set.K(data = data_1, LOCO = T, n.core = 60)
Yi_data_3 <- GWASpoly(data = data_2, models = models_1, traits = trait1, params = params, n.core = 50)

# save.image("~/Documents/git/big_files/data_Yi_DS_20.RData")
# load("~/Documents/git/big_files/data_Yi_DS_20.RData")


ggsave(filename = "~/Documents/git/big_files/M0_sqrt.jpg", plot = M0, width = 16, height = 16)

data_5 <- set.threshold(Yi_data_3, method= "Bonferroni", level=0.05)
QTL_01 <- get.QTL(data_5)
# manhattan.plot(data = data_5)


lev0 <- unique(QTL_01$Trait)

M0 <- manhattan.plot(data = data_5, traits = lev0) + theme_classic(base_family = "Arial", base_size = 12) + theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y = element_text(size = 12), plot.tag = element_blank()) 

ggsave(filename = "~/Documents/git/big_files/M0_sqrt_DS.jpg", plot = M0, width = 8, height = 8)

QTL_02 <- QTL_01 %>% distinct(QTL_01$Marker, .keep_all = T)



# Fit QTL -----------------------------------------------------------------

cc <- count(QTL_01,Trait)
lev4 <- cc$Trait
lev4
cc1 <- count(QTL_01, Model)
cc1$Model

QTL_03 <- QTL_01 %>% dplyr::filter(!Model %in% c("diplo-general", "diplo-additive"))


fit_05 <- fit.QTL(data=data_5, trait = "X20",
                  qtl=QTL_03[,c("Marker","Model")])


fit_06 <- list()
for (i in 1:length(lev4)) {
  fit_05 <- fit.QTL(data=data_5, trait = lev4[i],
                    qtl=QTL_03[,c("Marker","Model")])
  fit_05 <- fit_05 %>% group_by(Marker) %>% top_n(1, abs(R2))
  fit_06[[length(fit_06) + 1]] <- fit_05
}

names(fit_06) <- lev4
fit_06 <-rbindlist(fit_06, use.names=TRUE, fill=TRUE, idcol="trait")
fit_06 <- fit_06 %>% group_by(Marker) %>% top_n(1, abs(R2))
colnames(fit_06)
fit_06 <- fit_06[,c(2,6)]

fit_06$trait <- "Sum_Yi"
fit_06 <- fit_06 %>% distinct(Marker, .keep_all = T) 


# Tidy --------------------------------------------------------------------

QTL_03 <- QTL_01 %>% group_by(Marker) %>% top_n(1, abs(Score)) %>% dplyr::select(Marker, Score) %>% distinct(Marker, .keep_all = TRUE)
QTL_04 <- QTL_01 %>% group_by(Marker) %>% summarise(Trait = paste(Trait, collapse = ";")) 
QTL_05 <- QTL_01 %>% group_by(Marker) %>% summarise(Model = paste(Model, collapse = ";")) 
QTL_06 <- QTL_01 %>% dplyr::select(Marker, Chrom, Position) %>% distinct(Marker, .keep_all = TRUE)  %>% unite(col = "Marker1", 2:3, sep = "_", remove = T)

QTL_07 <-  QTL_06 %>% inner_join(., QTL_03, by = "Marker") %>% inner_join(., QTL_04, by = "Marker") %>% inner_join(., QTL_05, by = "Marker")

# annotate ----------------------------------------------------------------

load("~/Documents/git/big_files/i_5.2.8.RData")
anno2 <- read.table("~/Documents/git/big_files/annotations.txt", sep = "\t", header = T, fill = T, quote = "")
head(anno2)
dim(anno2)

file <- ("~/Documents/git/big_files/All_alfafa_pantr.bed")

txdb <- readBed(file, track.line = FALSE, remove.unusual = FALSE,
                zero.based = TRUE)

head(txdb)
mean(txdb$thickEnd-txdb$thickStart)


col_headings_2 <- c('gene_id',	'isoform')

QTL_02 <- QTL_01 %>% distinct(Marker, .keep_all = T) 
gr5 <- GRanges(seqnames = QTL_02$Chrom,
               ranges = IRanges(QTL_02$Position, width = 1))

overlaps <- join_overlap_left(gr5, txdb)

df2 <- annoGR2DF(overlaps)
df2 <- df2 %>% unite(col = "Marker1", 1:2, sep = "_", remove = T) %>% distinct(Marker1, .keep_all = TRUE) %>% dplyr::select(1,5) 
head(df2)

head(df3)
head(anno2)
df3 <- df2 %>% separate(2, col_headings_2, sep = ";", remove = TRUE, convert = FALSE, extra = "warn") %>% inner_join(., i_5.2.8, by = "gene_id")

df3 <- df2 %>% separate(2, col_headings_2, sep = ";", remove = TRUE, convert = FALSE, extra = "warn") %>% left_join(., anno2, by = "isoform")


QTL_08 <- inner_join(QTL_07, df3, by = "Marker1")

setwd("~/Documents/git/big_files/")
write.table(QTL_08, "markers_DS.tsv", row.names = F, quote = F, sep = "\t")


?GRanges
gr6 <- GRanges(seqnames = QTL_02$Chrom,
               ranges = IRanges(start=QTL_02$Position - 20000, end=QTL_02$Position + 20000, width=40001))


overlaps <- join_overlap_left(gr6, txdb)

df2 <- annoGR2DF(overlaps)
df2 <- df2 %>% unite(col = "Marker1", 1:2, sep = "_", remove = T) %>% distinct(name, .keep_all = TRUE) %>% dplyr::select(1,5) 
head(df2)

df3 <- df2 %>% separate(2, col_headings_2, sep = ";", remove = TRUE, convert = FALSE, extra = "warn") %>% left_join(., anno2, by = "isoform") %>% distinct(gene_id, .keep_all = TRUE)


setwd("~/Documents/git/big_files/")
write.csv(df3, "markers_DS1.csv", row.names = F, quote = F)
