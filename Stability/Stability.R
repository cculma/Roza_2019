rm(list = ls())

library(ASRtriala)
library(tidyverse)
library(svglite)
library(ggrepel)
setwd("~/Documents/git/big_files/")
C1 <- read.csv("BLUP_Yi_sqrt_SpATS_DArT.csv")
PCA <- read.csv("PCA_Roza2019.csv")
colnames(PCA)
PCA <- PCA %>% dplyr::select(Roza2019_VCF_ID, Plant_ID)

C1 <- inner_join(PCA, C1, by = "Roza2019_VCF_ID")
C1 <- C1[,-c(1,14,15,16)]
ncol(C1)
C2 <- C1 %>% gather(key = "Test", value = "predicted.value", 2:12)
C2$Test <- as.factor(C2$Test)
C2$Plant_ID <- as.factor(C2$Plant_ID)

summary(C2$Test)
C3 <- C2 %>% dplyr::filter(Test %in% c("yi_may", "yi_jun", "yi_jul", "yi_aug", "yi_sep"))
head(C3)
C3$g <- C3$predicted.value * 453.592
head(C3)
# Some stability tools
? ASRtriala::stability

# method = "wricke" (Wricke's ecovalence) 
# "static" (Shukla's stability variance)
# Shukla's stability variance (Shukla, 1972).
# "superiority" (cultivar-superiority measure),

stab.index <- stability(data = C3,
                        trial = "Test", gen = "Plant_ID", resp = "g",
                        method = "static", best = "max", plot = TRUE, top = TRUE,
                        bottom = FALSE, percentage = 5)
plot1 <- stab.index$stability.plot
plot1 <- plot1 + theme_classic(base_size = 14)
plot1
stab2 <- stab.index[["stats"]]
colnames(stab2)[2:3] <- c("static_ST2", "g_ST2")

setwd("~/medin297@umn.edu - Google Drive/My Drive/Roza2019/figs/")
svglite(filename = "Stability_ST2.svg", width = 5, height = 3, fix_text_size = F)
plot(plot1)
invisible(dev.off())


head(C3)
# Biplot
? ASRtriala::gbiplot
plot2 <- gbiplot(data = C3, vcov.g = NULL, scale = F,
        vector = "Test", unit = "Plant_ID", resp = "g",
        unit.label=TRUE)

plot2 <- plot2 + theme_classic(base_size = 14)
plot2


# ST1 ---------------------------------------------------------------------

setwd("~/Documents/git/big_files/")
D1 <- read.csv("BLUE_Yi_sqrt_SpATS_DArT.csv")
D1 <- inner_join(PCA, D1, by = "Roza2019_VCF_ID")
colnames(D1)
D1 <- D1[,-c(1,21:27)]
D2 <- D1 %>% gather(key = "Test", value = "predicted.value", 2:ncol(D1))
D2$Test <- as.factor(D2$Test)
D2$Plant_ID <- as.factor(D2$Plant_ID)
str(D2)
D2$g <- D2$predicted.value * 453.592
head(D2)

# "static" (Shukla's stability variance), 
stab.index <- stability(data = D2,
                        trial = "Test", gen = "Plant_ID", resp = "g",
                        method = "static", best = "max", plot = TRUE, top = TRUE,
                        bottom = FALSE, percentage = 5)
plot1 <- stab.index$stability.plot
plot1 <- plot1 + theme_classic(base_size = 14)
plot1

stab1 <- stab.index[["stats"]]
colnames(stab1)[2:3] <- c("static_ST1", "g_ST1")
stab3 <- inner_join(stab1,stab2, by = "Plant_ID")

setwd("~/Documents/git/Roza_2019/Stability/")
write.csv(stab3, "stab3.csv", quote = F, row.names = F)

setwd("~/medin297@umn.edu - Google Drive/My Drive/Roza2019/figs/")
svglite(filename = "Stability_ST1.svg", width = 5, height = 3, fix_text_size = F)
plot(plot1)
invisible(dev.off())

stab.index <- stability(data = D2,
                        trial = "Test", gen = "Plant_ID", resp = "predicted.value",
                        method = "rank", best = "max", plot = TRUE, top = TRUE,
                        bottom = FALSE, percentage = 10)

setwd("~/medin297@umn.edu - Google Drive/My Drive/Roza2019/GP1/")

marker <- read.csv("markers1.csv")
marker <- marker %>% separate(col = "Marker", into = c("chrom","pos"), sep = "_", remove = F)
class(marker$pos)
marker$pos <- as.numeric(marker$pos)
marker$pos1 <- marker$pos + 10000
marker <- marker %>% unite(col = "Marker1", c(2,4), sep = "_", remove = F)

marker$chrom <- as.factor(marker$chrom)
levels(marker$chrom)
marker$chrom <- recode_factor(marker$chrom,
                                chr1.1 = "1", chr2.1 = "2", chr3.1 = "3",
                                chr4.1 = "4", chr5.1 = "5", chr6.1 = "6",
                                chr7.1 = "7", chr8.1 = "8")

write.csv(marker, "markers2.csv", quote = F, row.names = F)


setwd("~/Documents/git/big_files/")

e1 <- read.csv("Effects_mrbean.csv")

head(e1)
str(e1)
e1$Trait <- as.factor(e1$Trait)
e1$acc_num <- as.factor(e1$acc_num)

# Biplot
? ASRtriala::gbiplot
plot2 <- gbiplot(data = e1, vcov.g = NULL, scale = T, weight = "weight",
                 vector = "Trait", unit = "acc_num", resp = "predicted.values",
                 unit.label=F)

plot2 <- plot2 + theme_classic(base_size = 14)
plot2

# subset low and high stability

head(D2)
D3 <- D2 %>% dplyr::filter(Plant_ID %in% c(515, 584, 423, 271)) %>% dplyr::filter(Test %in% c("may_20", "may_21", "jun_22", "may_23", "sep_20", "sep_21", "sep_22", "aug_23"))

D3 <- droplevels(D3)
levels(D3$Test)
D3$Test <- recode_factor(D3$Test,
                         "may_20" = 1, "may_21" = 1, "jun_22" = 1, "may_23" = 1, 
                         "sep_20" = 2, "sep_21" = 2, "sep_22" = 2, "aug_23" = 2)

DF <- data.frame(X = rep(1:4, each = 20),
                 Y = c(rnorm(20, 78, 1), rnorm(20, 76, 1),
                       rnorm(20, 80, 0.5), rnorm(20, 74, 2)))
Means <- DF %>% group_by(X) %>% summarize(Avg = mean(Y))

ggplot() + 
  geom_boxplot(data = D3, mapping = aes(x = Plant_ID, y = predicted.value, color = Plant_ID))


head(D3)
str(D3)
ggplot() + 
  geom_boxplot(data = D3, mapping = aes(x = Test, y = predicted.value, color = Plant_ID)) +
  geom_point(data = Means, mapping = aes(x = X, y = Avg)) +
  geom_line(data = Means, mapping = aes(x = X, y = Avg))


