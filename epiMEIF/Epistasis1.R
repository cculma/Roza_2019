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

geno2 <- inner_join(marker2, geno1, by = c("Chrom", "Position"))
geno2 <- geno2 %>% unite(col = "marker", 1:2, sep = "_", remove = T) %>% column_to_rownames("marker")

geno2 <- as.data.frame(t(geno2)) %>% rownames_to_column("gen")

geno3 <- inner_join(pheno, geno2, by = "gen")
colnames(geno3)[3:5] <- paste0("C",colnames(geno3)[3:5])
geno3 <- geno3 %>% dplyr::filter(gen != "UMN4351_38")

colnames(geno3)
str(geno3)
head(geno3)

# plot1 <- ggscatter(geno3, x = "C2_2963182", y = "winter_injury", size = 3, add = "reg.line", conf.int = T, point = T, alpha = 0.4) + stat_cor(label.x = 1, label.y = 5) + stat_regline_equation(label.x = 1, label.y = 4.5)
# plot1
# 
# plot2 <- ggscatter(geno3, x = "C6_11528532", y = "winter_injury", size = 3, add = "reg.line", conf.int = T, point = T, alpha = 0.4) + stat_cor(label.x = 1, label.y = 5) + stat_regline_equation(label.x = 1, label.y = 4.5)
# plot2
# 
# plot3 <- ggscatter(geno3, x = "C4_27397835", y = "winter_injury", size = 3, add = "reg.line", conf.int = T, point = T, alpha = 0.4) + stat_cor(label.x = 1, label.y = 5) + stat_regline_equation(label.x = 1, label.y = 4.5)
# plot3

geno4 <- geno3
geno4$C2_2963182 <- as.factor(geno4$C2_2963182)
geno4$C4_27397835 <- as.factor(geno4$C4_27397835)
geno4$C6_11528532 <- as.factor(geno4$C6_11528532)

ggplot(data = geno4, aes(x= C2_2963182, y=winter_injury)) + geom_boxplot() + theme_classic(base_size = 12) + theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.2))

ggplot(data = geno4, aes(x= C4_27397835, y=winter_injury)) + geom_boxplot() + theme_classic(base_size = 12) + theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.2))

ggplot(data = geno4, aes(x= C6_11528532, y=winter_injury)) + geom_boxplot() + theme_classic(base_size = 12) + theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.2))

head(geno4)

plot1 <- ggline(geno4, x = "C4_27397835", y = "winter_injury", add = "mean_se")+
  stat_compare_means() + 
  stat_compare_means(ref.group = "3", label = "p.signif", label.y = 4) + ylim(1, 5)  
plot1

plot2 <- ggline(geno4, x = "C6_11528532", y = "winter_injury", add = "mean_se")+
  stat_compare_means() +
  stat_compare_means(ref.group = "2", label = "p.signif", label.y = 4) + ylim(1, 5)    
plot2

# ggline(geno4, x = "C6_11528532", y = "winter_injury", add = "mean_se",linetype = "C4_27397835", shape = "C4_27397835")+
#   stat_compare_means() +
#   stat_compare_means(aes(group = C6_11528532), label = "p.signif", 
#                      label.y = c(4.5, 4.5, 4.5))


plot3 <- ggline(geno4, x = "C4_27397835", y = "winter_injury", add = "mean_se",linetype = "C6_11528532", shape = "C6_11528532")+
  stat_compare_means() +
  stat_compare_means(aes(group = C6_11528532), label = "p.signif", 
                     label.y = c(4.5, 4.5, 4.5)) + theme(legend.position = "none")  + ylim(1, 5)

plot3
cc <- dplyr::count(geno4, C4_27397835, C6_11528532)

setwd("~/medin297@umn.edu - Google Drive/My Drive/DAl21-6679/24_epistasis/")
write.csv(cc, "epistasis1.csv", quote = F, row.names = F)

plot4 <- ggarrange(plot1, plot2, plot3, ncol = 3, nrow = 1)

setwd("~/medin297@umn.edu - Google Drive/My Drive/DAl21-6679/4_figs/")
svglite(filename = "Epistasis1.svg", width = 6, height = 2.5, fix_text_size = F)
plot(plot4)
invisible(dev.off())

# ggplot(geno4, aes(x=C6_11528532, y=C4_27397835, size=winter_injury)) + geom_point()


plot5 <- ggline(geno4, x = "C4_27397835", y = "winter_injury", add = "mean_se",linetype = "C6_11528532", shape = "C6_11528532")+
  stat_compare_means() +
  stat_compare_means(aes(group = C6_11528532), label = "p.signif", 
                     label.y = c(4.5, 4.5, 4.5))  + ylim(1, 5)

setwd("~/medin297@umn.edu - Google Drive/My Drive/DAl21-6679/4_figs/")
svglite(filename = "Epistasis2.svg", width = 6, height = 2.5, fix_text_size = F)
plot(plot5)
invisible(dev.off())


plot5 <- ggline(geno4, x = "C4_27397835", y = "winter_injury", add = "mean_se",linetype = "C6_11528532", shape = "C6_11528532")+
  stat_compare_means() +
  stat_compare_means(aes(group = C6_11528532), label = "p.signif", 
                     label.y = c(4.5, 4.5, 4.5))  + ylim(1, 5)
