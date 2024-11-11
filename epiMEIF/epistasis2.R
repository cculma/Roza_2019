rm(list = ls())
library(GWASpoly)
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(svglite)

setwd("~/medin297@umn.edu - Google Drive/My Drive/DAl21-6679/3_GWASpoly/")
marker1 <- read.csv("GWAS1.csv")
marker2 <- marker1 %>% dplyr::filter(Trait == "Winter_injury")
marker2 <- marker2[c(3:5),c(5,6)]

pheno <- read.csv("Stem_strength_1.csv") # structure
head(pheno)
pheno <- pheno %>% dplyr::select(gen, winter_injury)

geno1 <- read.csv("DAl21_GWAS.csv", row.names = 1)
geno1[1:5,1:5]

geno2 <- geno1 %>% dplyr::filter(Chrom == 3) %>% dplyr::filter(Position > 77200000)
# top markers RF

marker3 <- data.frame(Chrom = c(4,2), Position = c(27397835, 2963182
))

# top markers SVM
marker4 <- data.frame(Chrom = c(4,5), Position = c(27397835, 35263239
))

geno2 <- inner_join(marker4, geno1, by = c("Chrom", "Position"))
geno2 <- geno2 %>% unite(col = "marker", 1:2, sep = "_", remove = T) %>% column_to_rownames("marker")

geno2 <- as.data.frame(t(geno2)) %>% rownames_to_column("gen")

geno3 <- inner_join(pheno, geno2, by = "gen")
colnames(geno3)
colnames(geno3)[3:ncol(geno3)] <- paste0("C",colnames(geno3)[3:ncol(geno3)])
geno3 <- geno3 %>% dplyr::filter(gen != "UMN4351_38")

geno4 <- geno3
head(geno4)
geno4$C4_27397835 <- as.factor(geno4$C4_27397835)
# geno4$C2_2963182 <- as.factor(geno4$C2_2963182)
geno4$C5_35263239 <- as.factor(geno4$C5_35263239)

# ggplot(data = geno4, aes(x= C2_2963182, y=winter_injury)) + geom_boxplot() + theme_classic(base_size = 12) + theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.2))

plot1 <- ggline(geno4, x = "C4_27397835", y = "winter_injury", add = "mean_se")+
  stat_compare_means() + 
  stat_compare_means(ref.group = "3", label = "p.signif", label.y = 4) + ylim(1, 5)  
plot1

plot2 <- ggline(geno4, x = "C5_35263239", y = "winter_injury", add = "mean_se")+
  stat_compare_means() +
  stat_compare_means(ref.group = "1", label = "p.signif", label.y = 4) + ylim(1, 5)    
plot2

# + theme(legend.position = "none")

plot3 <- ggline(geno4, x = "C4_27397835", y = "winter_injury", add = "mean_se",linetype = "C2_2963182", shape = "C2_2963182")+
  stat_compare_means() +
  stat_compare_means(aes(group = C2_2963182), label = "p.signif", 
                     label.y = c(4.5, 4.5, 4.5))  + ylim(1, 5)
plot3


plot3 <- ggline(geno4, x = "C4_27397835", y = "winter_injury", add = "mean_se",linetype = "C5_35263239", shape = "C5_35263239")+
  stat_compare_means() +
  stat_compare_means(aes(group = C5_35263239), label = "p.signif", 
                     label.y = c(4.5, 4.5, 4.5))  + ylim(1, 5)
plot3


# mixed linear model to epistasis


670.6 + 1785.6
959.6 + 1901.8
895.7 + 1516.1
1129.1 + 1541.6
753.6 + 1984.0
139.6 + 289.8 

1785.6 * 100 / (670.6 + 1785.6)
1901.8 * 100 / (959.6 + 1901.8)
1516.1 * 100 / (895.7 + 1516.1)
1541.6 * 100 / (1129.1 + 1541.6)
1984.0 * 100 / (753.6 + 1984.0)
289.8 * 100 / (139.6 + 289.8)
