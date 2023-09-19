# BLUP values of yield Roza2019
rm(list = ls())

library(SpATS)
library(tidyverse)
library(data.table)
library(rstatix)
library(stringr)

# library(ggcorrplot)
# library(asreml)
# library(asremlPlus)
# library(patchwork)

############

a1 <- read.csv("~/Documents/git/Roza_2019/raw_data/Roza2019_yield_raw.csv")
head(a1)
colnames(a1)
a1.1 <- a1[,c(2,4:8)] # position


a2 <- a1[,c(2,9:27)] # yield
a2 <- gather(a2, key = "harv", value = "yield", 2:20)
colnames(a2)
boxplot(a2$yield)
head(a2)

df_outliers <- a2 %>% group_by(harv) %>% identify_outliers(yield)
df_outliers1 <- df_outliers %>% unite(col = "id", 1:2, sep = "_", remove = T)
a3 <- a2 %>% unite(col = "id", 2:1, sep = "_", remove = F)
a4 <- a3 %>% dplyr::filter(!id %in% df_outliers1$id)
a5 <- a3 %>% dplyr::filter(id %in% df_outliers1$id)
a5$yield <- NA
a6 <- rbind(a4, a5)

boxplot(a6$yield)

a6 <- a6[,-1]

a7 <- inner_join(a1.1, a6, by = "plot")
head(a7)
ggplot(data = a7, aes(x = harv, y = yield)) + geom_boxplot() + theme_bw(base_size = 12, base_family = "Arial")

a7 <- inner_join(a1.1, a6, by = "plot")
a7 <- a7 %>% spread(key = harv, value = yield)
head(a7)

lev2 <- colnames(a7)[c(1,2,5,6)]
a7[,lev2] <- lapply(a7[,lev2], factor)

a7$acc_num <- as.factor(a7$acc_num)
a7$rep <- as.factor(a7$rep)
nlevels(a7$acc_num)
a7$R <- as.factor(a7$row)
a7$C <- as.factor(a7$col)
nlevels(a7$C)
nlevels(a7$R)
str(a7)
dim(a7)
lev3 <- colnames(a7)[7:25]
lev3

h2 <- list()
Y1 <- list()
for (i in 1:length(lev3)) {
  m1 <- SpATS(response = lev3[i],
              spatial = ~PSANOVA(col, row, nseg = c(nlevels(a7$C),nlevels(a7$R)), degree = c(3,3), nest.div = 2),
              fixed = ~ rep,
              random = ~ rep:C + rep:R,
              genotype.as.random = T,
              genotype = "acc_num",
              data = a7,
              control = list(tolerance = 1e-03))
  
  h1 <- getHeritability(m1) 
  pred4 <- predict.SpATS(m1, which = "acc_num", predFixed = "conditional")
  pred4 <- pred4[,c(1,7,8)]
  pred4$weight <- (1/pred4$standard.errors)^2
  
  h2[[length(h2)+1]] <- h1
  Y1[[length(Y1)+1]] <- pred4
}
names(Y1) <- lev3
names(h2) <- lev3
h3 <- as.data.frame(do.call(rbind, h2))
h3 <- h3 %>% rownames_to_column("env") 
colnames(h3)[2] <- "H2_SpATS"

setwd("~/Documents/git/Roza_2019/")
write.csv(h3, "H2_SpATS.csv", quote = F, row.names = F)


Y2 <-rbindlist(Y1, use.names=TRUE, fill=TRUE, idcol="env")
head(Y2)
colnames(Y2)[2] <- "gen"
Y2$env <- as.factor(Y2$env)
levels(Y2$env)

ggplot(data = Y2, aes(x = env, y = predicted.values)) + geom_boxplot() + theme_bw(base_size = 12, base_family = "Arial")

Y4 <- Y2 %>% dplyr::filter(!env %in% c("total_20","total_21","total_22"))
Y4$env <- droplevels(Y4$env)
levels(Y4$env)
head(Y4)
Y4 <- Y4 %>% separate(1, c("month", "year"), sep = "_", remove = F, convert = FALSE, extra = "merge")
lev5 <- c("month", "year")
Y4 <- as.data.frame(Y4)
Y4[,lev5] <- lapply(Y4[,lev5], factor)
Y4$month <- fct_relevel(Y4$month, c("may", "jun", "jul", "aug", "sep"))
str(Y4)

levels(Y4$month)
levels(Y4$year)
levels(Y4$gen)

data <- Y4
data <- data[order(data$gen, data$year), ]
data1 <- na.omit(data)
head(data1)
str(data1)

FA_1 <- asreml::asreml(fixed = predicted.values ~ 1, 
                       random = ~ sfa(year, 1):ar1(month):id(gen) + env,
                       data = data1, na.action = list(x = "include", y = "include"), 
                       weights = weight, family = asreml::asr_gaussian(dispersion = 1))


FA_1 <- update.asreml(FA_1)
summary(FA_1)$varcomp
current.asrt <- as.asrtests(FA_1, NULL, NULL)
current.asrt <- rmboundary.asrtests(current.asrt)


diffs <- predictPlus(classify = "gen:month", 
                     asreml.obj = FA_1, 
                     wald.tab = current.asrt$wald.tab, 
                     present = c("month","gen","year","env"))
ST2 <-diffs[[1]]

diffs <- predictPlus(classify = "gen:year", 
                     asreml.obj = FA_1, 
                     wald.tab = current.asrt$wald.tab, 
                     present = c("month","gen","year","env"))
ST3 <-diffs[[1]]

diffs <- predictPlus(classify = "gen", 
                     asreml.obj = FA_1, 
                     wald.tab = current.asrt$wald.tab, 
                     present = c("month","gen","year","env"))
ST4 <-diffs[[1]]

ST1 <- predict.asreml(FA_1, classify='gen:env', vcov=F)$pvals


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ST2$month <- paste0("ST2_", ST2$month)
ST3$year <- paste0("ST3_", ST3$year)

colnames(ST4)[2] <- "ST4_Yi"
ST4 <- ST4 %>% select(1:2)

ST2 <- ST2 %>% select(1:3) %>% spread(key = month, value = predicted.value) 
ST3 <- ST3 %>% select(1:3) %>% spread(key = year, value = predicted.value)
ST2 <- ST2[,c(1,5,4,3,2,6)]

head(Y2)
Y2.1 <- Y2
Y2.1$env <- paste0("ST0_", Y2.1$env)
Y2.1$env <- as.factor(Y2.1$env)
levels(Y2.1$env)
ST0 <- Y2.1 %>% select(1:3) %>% spread(key = env, value = predicted.values)


head(ST1)
ST1$env <- paste0("ST1_", ST1$env)
ST1 <- ST1 %>% select(1:3) %>% spread(key = env, value = predicted.value)

BLUP5 <- inner_join(ST0, ST1, by = "gen") %>% inner_join(., ST2, by = "gen") %>% inner_join(., ST3, by = "gen")  %>% inner_join(., ST4, by = "gen")
colnames(BLUP5)

write.csv(BLUP5, "~/Documents/git/Roza_2019/pheno_data/BLUEs_Yi_Roza2019.csv", quote = F, row.names = F)

#####################
# ST1 FA1
head(Y2)

Y2$env <- as.factor(Y2$env)
levels(Y2$env)


# Y4 <- Y2 %>% dplyr::filter(env %in% c("total_20","total_21","total_22"))

Y4 <- Y2 %>% dplyr::filter(!env %in% c("total_20","total_21","total_22"))
Y4$env <- droplevels(Y4$env)
levels(Y4$env)
head(Y4)
Y4 <- Y4 %>% separate(1, c("month", "year"), sep = "_", remove = F, convert = FALSE, extra = "merge")
lev5 <- c("month", "year")
Y4 <- as.data.frame(Y4)
Y4[,lev5] <- lapply(Y4[,lev5], factor)
Y4$month <- fct_relevel(Y4$month, c("may", "jun", "jul", "aug", "sep"))
str(Y4)
Y4$gen <- as.factor(Y4$gen)
levels(Y4$month)
levels(Y4$year)

data <- Y4
data <- data[order(data$gen, data$env), ]
data1 <- na.omit(data)
head(data1)
str(data1)



FA_1 <- update.asreml(FA_1)
summary(FA_1)$varcomp
current.asrt <- as.asrtests(FA_1, NULL, NULL)
current.asrt <- rmboundary.asrtests(current.asrt)

# asreml.options(workspace="128mb")
# asreml.options(workspace="512mb")
# asreml.options(workspace="1024mb")

# ~~~~~~~~~~~~~~~~~~~~~

# Y2 <- read.csv("~/Documents/git/Roza_2019/raw_data/BLUE_roza2019.csv")
# colnames(Y2)[1:2] <- c("env","gen") 
head(Y2)


Y3 <- Y2 %>% arrange(env) %>% mutate(env = factor(env, levels= lev3)) %>% select(1:3) %>% spread(key = env, value = predicted.values) %>% column_to_rownames("gen")

Y3 <- ST1 %>% unite("env", c(year, month), sep = "_", remove = F) %>% select(c(1,2,5)) %>% spread(key = env, value = predicted.value) %>% column_to_rownames("gen")


Y3 <- cor(Y3, use = "complete")

# ggplot(Y3, aes(x = env, y = predicted.values)) + geom_boxplot(outlier.shape = NA, alpha = 0.6, width=0.6, position = position_dodge(width=0.8, preserve = "single")) + theme_bw(base_family = "Arial") + theme(legend.position = "none", panel.spacing = unit(0.3, "lines"), strip.text.x = element_text(size = 10), axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.2), axis.title = element_text(size = 12)) + labs(title = "Roza2019 yield ST10", y = "", x = "") + ylim(0, 400)


save.image("~/Documents/git/big_files/SpATS1.RData")
load("~/Documents/git/big_files/SpATS1.RData")

# plots


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

