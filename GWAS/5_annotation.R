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

setwd("~/Documents/git/big_files/")
load("Yi_st2_51081_sqrt.RData")
load("Yi_st1_51081_sqrt.RData")
load("data_Yi_DS_20.RData")

data_5 <- set.threshold(Yi_data_3, method= "Bonferroni", level=0.05)
QTL_01 <- get.QTL(data_5)

# Tidy --------------------------------------------------------------------

QTL_03 <- QTL_01 %>% group_by(Marker) %>% top_n(1, abs(Score)) %>% dplyr::select(Marker, Score) %>% distinct(Marker, .keep_all = TRUE)
QTL_04 <- QTL_01 %>% group_by(Marker) %>% summarise(Trait = paste(Trait, collapse = ";")) 
QTL_05 <- QTL_01 %>% group_by(Marker) %>% summarise(Model = paste(Model, collapse = ";")) 
QTL_06 <- QTL_01 %>% dplyr::select(Marker, Chrom, Position) %>% distinct(Marker, .keep_all = TRUE)  %>% unite(col = "Marker1", 2:3, sep = "_", remove = T)
rm(QTL_07)
QTL_07 <-  QTL_06 %>% inner_join(., QTL_03, by = "Marker") %>% inner_join(., QTL_04, by = "Marker") %>% inner_join(., QTL_05, by = "Marker")

setwd("~/Documents/git/big_files/")
# write.table(QTL_07, "st2_51081_sqrt.tsv", row.names = F, quote = F, sep = "\t")

write.table(QTL_07, "st1_51081_sqrt.tsv", row.names = F, quote = F, sep = "\t")
QTL_07_6 <- QTL_07 # DS
QTL_07_7 <- QTL_07_6[,c(2:4)]

QTL_07_1 <- QTL_07 # st1
QTL_07_3 <- QTL_07_1[,c(2:4)]
QTL_07_2 <- QTL_07 # st2
QTL_07_4 <- QTL_07_2[,c(2:4)]

QTL_07_5 <- full_join(QTL_07_4,QTL_07_3, by = "Marker1")
colnames(QTL_07_5) <- c("Marker1", "Score_st2", "Trait_st2", "Score_st1", "Trait_st1")

QTL_07_8 <- full_join(QTL_07_5, QTL_07_7, by = "Marker1")
colnames(QTL_07_8)[6:7] <- c("Score_DS","Trait_DS")

setwd("~/Documents/git/big_files/")
write.table(QTL_07_8, "markers_join.tsv", row.names = F, quote = F, sep = "\t")

head(QTL_07_8)

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

QTL_02 <- QTL_07_8 %>% separate(1, c("chrom", "pos"), sep = "_", remove = TRUE, convert = FALSE, extra = "warn") 
QTL_02$chrom <- as.factor(QTL_02$chrom)
QTL_02$pos <- as.numeric(QTL_02$pos)

# genes in the same marker position
# gr5 <- GRanges(seqnames = QTL_02$chrom,
#                ranges = IRanges(QTL_02$pos, width = 1))
# overlaps <- join_overlap_left(gr5, txdb)

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

# using a range of 10 kb upstream and 10 kb downstream
gr6 <- GRanges(seqnames = QTL_02$chrom,
               ranges = IRanges(start=QTL_02$pos - 10000, end=QTL_02$pos + 10000, width=20001))
gr6

overlaps <- join_overlap_left(gr6, txdb)

df2 <- annoGR2DF(overlaps)
df2 <- df2 %>% unite(col = "Marker1", 1:2, sep = "_", remove = T) %>% distinct(name, .keep_all = TRUE) %>% dplyr::select(1,5) 
head(df2)

df3 <- df2 %>% separate(2, col_headings_2, sep = ";", remove = TRUE, convert = FALSE, extra = "warn") %>% left_join(., anno2, by = "isoform") %>% distinct(gene_id, .keep_all = TRUE)
head(df3)

setwd("~/Documents/git/big_files/")
# write.csv(df3, "markers_DS1.csv", row.names = F, quote = F)
# write.csv(df3, "markers_Yi_GO.csv", row.names = F, quote = F)

GO1 <- read.csv("markers_Yi_GO1.csv")

lev4 <- c("GO1","GO2","GO3","GO4","GO5","GO6","GO7","GO8")
GO2 <- GO1 %>% separate(3, lev4, sep = ";", remove = TRUE, convert = FALSE, extra = "warn") 

GO2 <- GO2 %>% gather(key = "GO", value = "term", 3:10) 
GO2$term <- gsub("^\\s+|\\s+$","",GO2$term)


GO2 <- GO2 %>% separate(4, c("term","GO_ID"), sep = "\\[G", remove = TRUE, convert = FALSE, extra = "warn") 

GO2$GO_ID <- gsub("^[^_]*_|\\]","",GO2$GO_ID)

cc <- dplyr::count(GO2, GO_ID)
cc <- dplyr::count(GO2, term)

setwd("~/Documents/git/big_files/")
write.table(cc, "count_GO.tsv", row.names = F, quote = F, sep = "\t")
write.table(GO2, "markers_Yi_GO2.tsv", row.names = F, quote = F, sep = "\t")
