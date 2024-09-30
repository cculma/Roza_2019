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

load("~/Documents/git/big_files/Yi_st1_51081_sqrt.RData")
data_3 <- set.threshold(Yi_data_3, method= "Bonferroni", level=0.05)
QTL_03 <- get.QTL(data_3)

load("~/Documents/git/big_files/Yi_st2_51081_sqrt.RData")
data_4 <- set.threshold(Yi_data_3, method= "Bonferroni", level=0.05)
QTL_04 <- get.QTL(data_4)

load("~/Documents/git/big_files/data_Yi_DS_20.RData")
data_5 <- set.threshold(Yi_data_3, method= "Bonferroni", level=0.05)
QTL_05 <- get.QTL(data_5)

QTL_06 <- rbind(QTL_03, QTL_04, QTL_05)
QTL_06.1 <- QTL_06 %>% distinct(Marker, .keep_all = T)


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


QTL_06 <- rbind(QTL_03, QTL_04, QTL_05)
QTL_06.1 <- QTL_06 %>% distinct(Marker, .keep_all = T)

QTL_06.1$Chrom <- as.factor(QTL_06.1$Chrom)
QTL_06.1$Position <- as.numeric(QTL_06.1$Position)

# genes in the same marker position

gr5 <- GRanges(seqnames = QTL_06.1$Chrom,
               ranges = IRanges(QTL_06.1$Position, width = 1))
overlaps <- join_overlap_left(gr5, txdb)

df2 <- annoGR2DF(overlaps)
df2 <- df2 %>% unite(col = "Marker", 1:2, sep = "_", remove = T) %>% distinct(Marker, .keep_all = TRUE) %>% dplyr::select(1,5) 
head(df2)

head(i_5.2.8)
head(anno2)
df3 <- df2 %>% separate(2, col_headings_2, sep = ";", remove = TRUE, convert = FALSE, extra = "warn") %>% inner_join(., i_5.2.8, by = "gene_id")
head(df3)

df3 <- df2 %>% separate(2, col_headings_2, sep = ";", remove = TRUE, convert = FALSE, extra = "warn") %>% left_join(., anno2, by = "isoform")

df3 <- df2 %>% separate(2, col_headings_2, sep = ";", remove = TRUE, convert = FALSE, extra = "warn") %>% inner_join(., anno2, by = "isoform")

df4 <- na.omit(df3) 

QTL_08 <- inner_join(QTL_07, df3, by = "Marker1")

setwd("~/Documents/git/big_files/")
write.table(QTL_08, "markers_DS.tsv", row.names = F, quote = F, sep = "\t")

# using a range of 10 kb upstream and 10 kb downstream
# using a range of 1 kb upstream and 1 kb downstream
gr6 <- GRanges(seqnames = QTL_06.1$Chrom,
               ranges = IRanges(start=QTL_06.1$Position - 1000, end=QTL_06.1$Position + 1000, width=2001))
gr6

overlaps <- join_overlap_left(gr6, txdb)

df2 <- annoGR2DF(overlaps)

df2 <- df2 %>% separate(col = "name", into = c("gene", "isoform"), sep = ";", remove = T) %>% distinct(gene, .keep_all = TRUE)

df2$Position <- df2$start + 1000

df3 <- df2 %>% left_join(., anno2, by = "isoform")
colnames(df3)
head(df3)

# df3 <- df3 %>% dplyr::select("chr","thickStart","thickEnd","Position","Uniprot")
# 
# colnames(df3)[1] <- "Chrom"
df3 <- df3 %>% unite(col = "Marker", sep = "_", c("chr","Position"), remove = F)

# df3 <- df2 %>% separate(2, col_headings_2, sep = ";", remove = TRUE, convert = FALSE, extra = "warn") %>% left_join(., anno2, by = "isoform") %>% distinct(gene_id, .keep_all = TRUE)
head(df3)
dim(df3)

df4 <- df3 %>% dplyr::select(Marker, Uniprot)
colnames(df4)[2] <- "UniProt"
df4 <- na.omit(df4)
head(df4)
df4$new <- 1

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

setwd("~/Documents/git/big_files/")
anno1 <- read.csv("roza_2019_anno.csv")
head(anno1)
dim(anno1)

head(df3)
dim(df3)
colnames(df3)[1] <- "Marker"

anno3 <- inner_join(anno1, df4, by = c("Marker","UniProt"))

anno3 <- left_join(anno1, df4, by = c("Marker","UniProt"))
setwd("~/Documents/git/big_files/")
write.table(anno3, "roza_2019_anno1.csv", row.names = F, quote = F, sep = "\t")
