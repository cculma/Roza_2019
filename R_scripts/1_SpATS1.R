# BLUP values of yield Roza2019
rm(list = ls())

library(SpATS)
library(tidyverse)
library(data.table)
library(ggcorrplot)
library(asreml)
library(asremlPlus)
library(patchwork)

############


a1 <- read.csv("~/Documents/git/Roza_2019/raw_data/Roza2019_yield.csv")
head(a1)
colnames(a1)
lev2 <- c("acc","acc_num","rep")
a1[,lev2] <- lapply(a1[,lev2], factor)

nlevels(a1$acc)
nlevels(a1$acc_num)
a1$R <- as.factor(a1$row)
a1$C <- as.factor(a1$col)
nlevels(a1$C)
nlevels(a1$R)

lev3 <- colnames(a1)[7:23]

lev4 <- c("may_20","jun_20","jul_20","aug_20","sep_20","total_20",
          "may_21","jun_21","jul_21","aug_21","sep_21","total_21",
          "jun_22","jul_22","aug_22","sep_22","total_22")


Y1 <- list()
for (i in 1:length(lev3)) {
  m1 <- SpATS(response = lev3[i],
              spatial = ~PSANOVA(col, row, nseg = c(nlevels(a1$C),nlevels(a1$R)), degree = c(3,3), nest.div = 2),
              fixed = ~rep,
              random = ~ rep:C + rep:R,
              genotype.as.random = F,
              genotype = "acc_num",
              data = a1,
              control = list(tolerance = 1e-03))
  
  pred4 <- predict.SpATS(m1, which = "acc_num", predFixed = "marginal")
  pred4 <- pred4[,c(1,7,8)]
  pred4$weight <- (1/pred4$standard.errors)^2
  Y1[[length(Y1)+1]] <- pred4
}
names(Y1) <- lev3

Y2 <-rbindlist(Y1, use.names=TRUE, fill=TRUE, idcol="env")
head(Y2)

Y3 <- Y2 %>% arrange(env) %>% mutate(env = factor(env, levels= lev4)) %>% select(1:3) %>% spread(key = env, value = predicted.values) %>% column_to_rownames("acc_num")
Y3 <- cor(Y3, use = "complete")

# ggplot(Y3, aes(x = env, y = predicted.values)) + geom_boxplot(outlier.shape = NA, alpha = 0.6, width=0.6, position = position_dodge(width=0.8, preserve = "single")) + theme_bw(base_family = "Arial") + theme(legend.position = "none", panel.spacing = unit(0.3, "lines"), strip.text.x = element_text(size = 10), axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.2), axis.title = element_text(size = 12)) + labs(title = "Roza2019 yield ST10", y = "", x = "") + ylim(0, 400)


# save.image("~/Documents/git/big_files/SpATS1.RData")
# load("~/Documents/git/big_files/SpATS1.RData")
#####################
# ST1 FA1
head(Y2)
colnames(Y2)[2] <- "gen" 
Y2$env <- as.factor(Y2$env)
levels(Y2$env)

Y4 <- Y2 %>% dplyr::filter(!env %in% c("total_20","total_21","total_22"))
# Y4 <- Y2 %>% dplyr::filter(env %in% c("total_20","total_21","total_22"))
Y4$env <- droplevels(Y4$env)
levels(Y4$env)
head(Y4)
Y4 <- Y4 %>% separate(1, c("month", "year"), sep = "_", remove = F, convert = FALSE, extra = "merge")
str(Y4)
lev5 <- c("month", "year")
Y4 <- as.data.frame(Y4)
Y4[,lev5] <- lapply(Y4[,lev5], factor)

# data <- Y2
data <- Y2
data <- data[order(data$gen, data$env), ]
data1 <- na.omit(data)
head(data1)
str(data1)

FA_1 <- asreml::asreml(fixed = predicted.values ~ 1, 
                       random = ~ fa(env, 2):id(gen),
                       data = data1, na.action = list(x = "include", y = "include"), 
                       weights = weight, family = asreml::asr_gaussian(dispersion = 1))

FA_1 <- update.asreml(FA_1)


FA_1 <- asreml::asreml(fixed = predicted.values ~ 1 + gen, 
                       random = ~ fa(env, 2):id(gen) + month + year,
                       data = data1, na.action = list(x = "include", y = "include"), 
                       weights = weight, family = asreml::asr_gaussian(dispersion = 1))


# asreml.options(workspace="128mb")
# asreml.options(workspace="512mb")
# asreml.options(workspace="1024mb")

# save.image("~/Documents/git/big_files/1_SpATS1.RData")
# load("~/Documents/git/big_files/1_SpATS1.RData")

ST1 <- predict.asreml(FA_1, classify='gen:env', vcov=F)$pvals
ST2 <- predict.asreml(FA_1, classify='gen:month', vcov=F)$pvals
ST3 <- predict.asreml(FA_1, classify='gen:year', vcov=F)$pvals
ST4 <- predict.asreml(FA_1, classify='gen', vcov=F)$pvals

hist(ST1$predicted.value)
hist(ST2$predicted.value)
hist(ST3$predicted.value)
hist(ST4$predicted.value)

head(BLUP2)

levels(BLUP2$env)
BLUP2$env <- factor(BLUP2$env, levels = lev4)

# BLUP2 <- BLUP2 %>% arrange(env) %>% mutate(env = factor(env, levels= lev4))

P2 <- ggplot(BLUP2, aes(x = env, y = predicted.value)) + geom_boxplot(outlier.shape = NA, alpha = 0.6, width=0.6, position = position_dodge(width=0.8, preserve = "single")) + theme_bw(base_family = "Arial", base_size = 12) + theme(legend.position = "none", panel.spacing = unit(0.3, "lines"), axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.2), axis.title = element_text(size = 12)) + labs(title = "Roza2019 yield ST1", y = "BLUP Yield (lb * plot * 100)", x = "") + ylim(0, 400)


head(Y2)
levels(Y2$env)
Y2$env <- factor(Y2$env, levels = lev4)

P3 <- ggplot(Y2, aes(x = env, y = predicted.values)) + geom_boxplot(outlier.shape = NA, alpha = 0.6, width=0.6, position = position_dodge(width=0.8, preserve = "single")) + theme_bw(base_family = "Arial", base_size = 12) + theme(legend.position = "none", panel.spacing = unit(0.3, "lines"), axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.2), axis.title = element_text(size = 12)) + labs(title = "Roza2019 yield ST0", y = "BLUE Yield (lb * plot * 100)", x = "") + ylim(0, 400)

BLUP4 <- BLUP2 %>% select(1:3) %>% spread(key = env, value = predicted.value) %>% column_to_rownames("gen")
BLUP4 <- cor(BLUP4, use = "complete")


Y3[lower.tri(Y3)] <- BLUP4[lower.tri(BLUP4)]

P1 <- ggcorrplot(Y3[,ncol(Y3):1], hc.order = F, type = "full", lab = T, lab_col = "grey3", lab_size = 3, show.diag = T) + theme_classic(base_family = "Arial", base_size = 12) + theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.2), axis.title.x=element_blank(), axis.title.y = element_blank()) + labs(title = "Single-Stage vs Stage-Wise")

# P2 + P1 + plot_layout(ncol = 2)
# 
# patch <- p1 + p2
# 
# patch <- P3 / P2
# patch + P1 + plot_layout(ncol = 2)
# 
# (P3 / P2) + P1 + plot_layout(ncol = 2)
# P3 / (P1 | P2)

# 12 x 8

