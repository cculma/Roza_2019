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
library(grid)
library(gtable)
library(plyranges)
library(Repitools)
library(devtools)
library(svglite)

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

QTL_05.20 <- QTL_05 %>% dplyr::filter(Trait == "X20")
QTL_05.21 <- QTL_05 %>% dplyr::filter(Trait == "X21")
QTL_05.22 <- QTL_05 %>% dplyr::filter(Trait == "X22")
QTL_05.23 <- QTL_05 %>% dplyr::filter(Trait == "X23")

manhattan.plot(data_5, traits = "X21",chrom = "chr2.1")
manhattan.plot(data_5, traits = "X22",chrom = "chr2.1")


dplyr::count(QTL_05.20, Model)
dplyr::count(QTL_05.21, Model)
dplyr::count(QTL_05.22, Model)
dplyr::count(QTL_05.23, Model)

QTL_06 <- rbind(QTL_03, QTL_04, QTL_05)
QTL_06.1 <- QTL_06 %>% distinct(Marker, .keep_all = T)
QTL_06.1$Model

cols <- c("1" = "#377EB8", "2" = "#636363", "3" = "#377EB8", "4" = "#636363",
          "5" = "#377EB8", "6" = "#636363", "7" = "#377EB8", "8" = "#636363")

scores <- data_5@scores[["X20"]]
scores <- data_5@scores[["X21"]]
scores <- data_5@scores[["X22"]]
scores <- data_5@scores[["X23"]]

rownames(scores) <- paste0(data_5@map$Chrom,  "_",data_5@map$Position)

scores <- scores %>% dplyr::select(`general`, `1-dom-alt`, `2-dom-ref`)
scores <- scores %>% dplyr::select(`1-dom-alt`, `1-dom-ref`, `2-dom-alt`, `diplo-general`, `general`)
scores <- scores %>% dplyr::select(`1-dom-alt`, `1-dom-ref`, `2-dom-ref`, `diplo-general`, `general`)
scores <- scores %>% dplyr::select(`1-dom-alt`, `1-dom-ref`, `2-dom-alt`, `diplo-general`, general)


G1.1 <- scores %>% rownames_to_column("Marker") %>% separate(col = "Marker", into = c("chrom", "pos"), sep = "_", remove = F)
G1.1$chrom <- as.factor(G1.1$chrom)
G1.1$Marker <- as.factor(G1.1$Marker)
G1.1$pos <- as.numeric(G1.1$pos)
levels(G1.1$chrom)
G1.1$chrom <- recode_factor(G1.1$chrom,
                            chr1.1 = "1", chr2.1 = "2", chr3.1 = "3",
                            chr4.1 = "4", chr5.1 = "5", chr6.1 = "6",
                            chr7.1 = "7", chr8.1 = "8")

G1.2 <- G1.1 %>% gather(key = model, value = p_val, 4:ncol(G1.1))
head(G1.2)
G1.2 <- na.omit(G1.2)
G1.2 <- G1.2 %>% dplyr::filter(p_val > 2)

plot1.2 <- ggplot(data = G1.2, aes(x = Marker, y = p_val, group = 1)) + geom_point(aes(colour = factor(chrom), alpha = 0.5), size = 1) + theme_classic(base_family = "Arial", base_size = 12) + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") + scale_color_manual(values = cols) + ylim(2, 9) + geom_hline(yintercept=5.8, linetype="dashed", color = "gray") + labs(y = "2020")
plot1.2

plot1.3 <- ggplot(data = G1.2, aes(x = Marker, y = p_val, group = 1)) + geom_point(aes(colour = factor(chrom), alpha = 0.5), size = 1) + theme_classic(base_family = "Arial", base_size = 12) + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") + scale_color_manual(values = cols) + ylim(2, 9) + geom_hline(yintercept=5.8, linetype="dashed", color = "gray") + labs(y = "2021")

plot1.4 <- ggplot(data = G1.2, aes(x = Marker, y = p_val, group = 1)) + geom_point(aes(colour = factor(chrom), alpha = 0.5), size = 1) + theme_classic(base_family = "Arial", base_size = 12) + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") + scale_color_manual(values = cols) + ylim(2, 9) + geom_hline(yintercept=5.8, linetype="dashed", color = "gray") + labs(y = "2022")

plot1.5 <- ggplot(data = G1.2, aes(x = Marker, y = p_val, group = 1)) + geom_point(aes(colour = factor(chrom), alpha = 0.5), size = 1) + theme_classic(base_family = "Arial", base_size = 12) + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") + scale_color_manual(values = cols) + ylim(2, 9) + geom_hline(yintercept=5.8, linetype="dashed", color = "gray") + labs(y = "2023")

plot1.2
plot1.3
plot1.4
plot1.5

plot1.6 <- ggarrange(plot1.2, plot1.3, plot1.4, plot1.5, ncol = 1, nrow = 4)

setwd("~/medin297@umn.edu - Google Drive/My Drive/Roza2019/figs/")
svglite(filename = "Manhattan_DS1.svg", width = 5, height = 4, fix_text_size = F)
plot(plot1.6)
invisible(dev.off())


# Figure S3 ---------------------------------------------------------------

QTL_03$Trait

trait1 <- QTL_03 %>% distinct(Trait, .keep_all = F)
lev1 <- trait1$Trait
trait2 <- QTL_04 %>% distinct(Trait, .keep_all = F)
trait2$Trait <- as.factor(trait2$Trait)
levels(trait2$Trait) 

trait2 <- trait2 %>% dplyr::filter(!Trait %in% c("yi_year", "yi_2020", "yi_2021", "yi_2022", "yi_2023"))
lev2 <- trait2$Trait
lev2 <- droplevels(lev2)

i = 3
lev3 <- list()
for (i in 1:length(lev1)) {
  
  QTL_05.20 <- QTL_03 %>% dplyr::filter(Trait == lev1[[i]])
  models_2 <- QTL_05.20 %>% distinct(Model, .keep_all = F)

  scores <- data_3@scores[[lev1[[i]]]]
  scores <- scores %>% dplyr::select(models_2$Model)
  rownames(scores) <- paste0(data_3@map$Chrom,  "_",data_3@map$Position)
  
  G1.1 <- scores %>% rownames_to_column("Marker") %>% separate(col = "Marker", into = c("chrom", "pos"), sep = "_", remove = F)
  G1.1$chrom <- as.factor(G1.1$chrom)
  G1.1$Marker <- as.factor(G1.1$Marker)
  G1.1$pos <- as.numeric(G1.1$pos)
  levels(G1.1$chrom)
  G1.1$chrom <- recode_factor(G1.1$chrom,
                              chr1.1 = "1", chr2.1 = "2", chr3.1 = "3",
                              chr4.1 = "4", chr5.1 = "5", chr6.1 = "6",
                              chr7.1 = "7", chr8.1 = "8")
  
  G1.2 <- G1.1 %>% gather(key = model, value = p_val, 4:ncol(G1.1))
  head(G1.2)
  G1.2 <- na.omit(G1.2)
  
  plot1.2 <- ggplot(data = G1.2, aes(x = Marker, y = p_val, group = 1)) + geom_point(aes(colour = factor(chrom), alpha = 0.5), size = 1) + theme_classic(base_family = "Arial", base_size = 12) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") + scale_color_manual(values = cols) + ylim(2, max(G1.2$p_val + 1)) + geom_hline(yintercept=5.8, linetype="dashed", color = "gray") + labs(x = lev1[[i]], y = expression(-log[10](p)))
  
  lev3[[length(lev3)+1]] <- plot1.2
}

names(lev3) <- lev1


lev4 <- list()
for (i in 1:length(lev2)) {
  
  QTL_05.20 <- QTL_04 %>% dplyr::filter(Trait == lev2[[i]])
  models_2 <- QTL_05.20 %>% distinct(Model, .keep_all = F)
  
  scores <- data_4@scores[[lev2[[i]]]]
  scores <- scores %>% dplyr::select(models_2$Model)
  rownames(scores) <- paste0(data_4@map$Chrom,  "_",data_4@map$Position)
  
  G1.1 <- scores %>% rownames_to_column("Marker") %>% separate(col = "Marker", into = c("chrom", "pos"), sep = "_", remove = F)
  G1.1$chrom <- as.factor(G1.1$chrom)
  G1.1$Marker <- as.factor(G1.1$Marker)
  G1.1$pos <- as.numeric(G1.1$pos)
  levels(G1.1$chrom)
  G1.1$chrom <- recode_factor(G1.1$chrom,
                              chr1.1 = "1", chr2.1 = "2", chr3.1 = "3",
                              chr4.1 = "4", chr5.1 = "5", chr6.1 = "6",
                              chr7.1 = "7", chr8.1 = "8")
  
  G1.2 <- G1.1 %>% gather(key = model, value = p_val, 4:ncol(G1.1))
  # head(G1.2)
  G1.2 <- na.omit(G1.2)
  
  plot1.2 <- ggplot(data = G1.2, aes(x = Marker, y = p_val, group = 1)) + geom_point(aes(colour = factor(chrom), alpha = 0.5), size = 1) + theme_classic(base_family = "Arial", base_size = 12) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") + scale_color_manual(values = cols) + ylim(2, max(G1.2$p_val + 1)) + geom_hline(yintercept=5.8, linetype="dashed", color = "gray") + labs(x = lev2[[i]], y = expression(-log[10](p)))
  
  lev4[[length(lev4)+1]] <- plot1.2
}

names(lev4) <- lev2
lev4[[1]]


lev5 <- c(lev3, lev4)
length(lev5)
names(lev5)


g1 <- ggplotGrob(lev5[[1]])
g2 <- ggplotGrob(lev5[[2]])
g3 <- ggplotGrob(lev5[[3]])
g4 <- ggplotGrob(lev5[[4]])
g5 <- ggplotGrob(lev5[[5]])

g6 <- ggplotGrob(lev5[[6]])
g7 <- ggplotGrob(lev5[[7]])
g8 <- ggplotGrob(lev5[[8]])
g9 <- ggplotGrob(lev5[[9]])
g10 <- ggplotGrob(lev5[[10]])

g11 <- ggplotGrob(lev5[[11]])
g12 <- ggplotGrob(lev5[[12]])
g13 <- ggplotGrob(lev5[[13]])
g14 <- ggplotGrob(lev5[[14]])
g15 <- ggplotGrob(lev5[[15]])
g16 <- ggplotGrob(lev5[[16]])

g17 <- ggplotGrob(lev5[[21]])
g18 <- ggplotGrob(lev5[[22]])
g19 <- ggplotGrob(lev5[[23]])
g20 <- ggplotGrob(lev5[[24]])
g21 <- ggplotGrob(lev5[[25]])
g22 <- ggplotGrob(lev5[[26]])

# 22 plots 6 rows 4 cols

c1 <- cbind(g1, g2, g3, g4, g5, g6, size = "first")
c2 <- cbind(g7, g8, g9, g10, g11, g12, size = "first")
c3 <- cbind(g13, g14, g15, g16, g17, g18, size = "first")
c4 <- cbind(g19, g20, g21, g22, g22, g22, size = "first")

c1 <- cbind(g1, g2, g3, g4, size = "first")
c2 <- cbind(g5, g6, g7, g8, size = "first")
c3 <- cbind(g9, g10, g11, g12, size = "first")
c4 <- cbind(g13, g14, g15, g16,  size = "first")
c5 <- cbind(g17, g18, g19, g20,  size = "first")
c6 <- cbind(g21, g22, g22, g22,  size = "first")

plot1.3 <- rbind(c1, c2, c3, c4, c5, c6, size = "first")

plot1.3 <- ggarrange(
          lev5[[1]], lev5[[2]], lev5[[3]], lev5[[4]], lev5[[5]],
          lev5[[6]], lev5[[7]], lev5[[8]], lev5[[9]], lev5[[10]],
          lev5[[11]], lev5[[12]], lev5[[13]], lev5[[14]], lev5[[15]],
          lev5[[16]], 
          lev5[[21]], lev5[[22]], lev5[[23]], lev5[[24]], lev5[[25]],
          lev5[[26]],
          ncol = 3, nrow = 9)


setwd("~/medin297@umn.edu - Google Drive/My Drive/Roza2019/figs/")
svglite(filename = "Manhattan_ST1.svg", width = 6, height = 8, fix_text_size = F)

plot(plot1.3)
invisible(dev.off())

