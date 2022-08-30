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

a1$acc <- as.factor(a1$acc)
a1$acc_num <- as.factor(a1$acc_num)
nlevels(a1$acc)
nlevels(a1$acc_num)
a1$R <- as.factor(a1$row)
a1$C <- as.factor(a1$col)
a1$rep <- as.factor(a1$rep)
nlevels(a1$C)
nlevels(a1$R)

lev3 <- colnames(a1)[7:22]

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

Y3[lower.tri(Y3)] <- BLUP4[lower.tri(BLUP4)]

P1 <- ggcorrplot(Y3[,ncol(Y3):1], hc.order = F, type = "full", lab = T, lab_col = "grey3", lab_size = 3, show.diag = T) + theme_classic(base_family = "Arial", base_size = 12) + theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.2), axis.title.x=element_blank(), axis.title.y = element_blank()) + labs(title = "Single-Stage vs Stage-Wise")

# ggplot(Y3, aes(x = env, y = predicted.values)) + geom_boxplot(outlier.shape = NA, alpha = 0.6, width=0.6, position = position_dodge(width=0.8, preserve = "single")) + theme_bw(base_family = "Arial") + theme(legend.position = "none", panel.spacing = unit(0.3, "lines"), strip.text.x = element_text(size = 10), axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.2), axis.title = element_text(size = 12)) + labs(title = "Roza2019 yield ST10", y = "", x = "") + ylim(0, 400)


# save.image("~/Documents/git/big_files/SpATS1.RData")
# load("~/Documents/git/big_files/SpATS1.RData")
#####################
# ST1 FA1
head(Y2)
colnames(Y2)[2] <- "gen" 
Y2$env <- as.factor(Y2$env)
levels(Y2$env)

# Y4 <- Y2 %>% dplyr::filter(!env %in% c("all_20","all_21","all_22"))
# Y4 <- Y2 %>% dplyr::filter(env %in% c("all_20","all_21","all_22"))
# Y4$env <- droplevels(Y4$env)
# levels(Y4$env)

data <- Y2
data <- data[order(data$gen, data$env), ]
data1 <- na.omit(data)
head(data1)
str(data1)

FA_1 <- asreml::asreml(fixed = predicted.values ~ 1, 
                       random = ~ fa(env, 1):id(gen),
                       data = data1, na.action = list(x = "include", y = "include"), 
                       weights = weight, family = asreml::asr_gaussian(dispersion = 1))

FA_1 <- update.asreml(FA_1)

# asreml.options(workspace="128mb")
# asreml.options(workspace="512mb")
# asreml.options(workspace="1024mb")

BLUP2 <- predict.asreml(FA_1, classify='gen:env', vcov=F)$pvals

head(BLUP2)

levels(BLUP2$env)
lev4 <- c("may_20","jun_20","jul_20","aug_20","sep_20","all_20",
          "may_21","jun_21","jul_21","aug_21","sep_21","all_21",
          "jun_22","jul_22","aug_22","all_22")

BLUP2 <- BLUP2 %>% arrange(env) %>% mutate(env = factor(env, levels= lev4))

P2 <- ggplot(BLUP2, aes(x = env, y = predicted.value)) + geom_boxplot(outlier.shape = NA, alpha = 0.6, width=0.6, position = position_dodge(width=0.8, preserve = "single")) + theme_bw(base_family = "Arial", base_size = 12) + theme(legend.position = "none", panel.spacing = unit(0.3, "lines"), axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.2), axis.title = element_text(size = 12)) + labs(title = "Roza2019 yield ST1", y = "", x = "") + ylim(0, 400)


BLUP4 <- BLUP2 %>% select(1:3) %>% spread(key = env, value = predicted.value) %>% column_to_rownames("gen")
BLUP4 <- cor(BLUP4, use = "complete")

P2 + P1 + plot_layout(ncol = 2)






