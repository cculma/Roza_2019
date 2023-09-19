# H2 values of yield Roza2019
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
a1.1 <- a1[,c(2,4:8)] # position info


a2 <- a1[,c(2,9:27)] # yield data
a2 <- gather(a2, key = "harv", value = "yield", 2:20)
colnames(a2)
boxplot(a2$yield)
head(a2)

# remove outliers

df_outliers <- a2 %>% group_by(harv) %>% identify_outliers(yield)
df_outliers1 <- df_outliers %>% unite(col = "id", 1:2, sep = "_", remove = T)
a3 <- a2 %>% unite(col = "id", 2:1, sep = "_", remove = F)
a4 <- a3 %>% dplyr::filter(!id %in% df_outliers1$id)
a5 <- a3 %>% dplyr::filter(id %in% df_outliers1$id)
a5$yield <- NA
a6 <- rbind(a4, a5)

boxplot(a6$yield)
a6 <- a6[,-1]

# a7 <- inner_join(a1.1, a6, by = "plot")
# head(a7)
# ggplot(data = a7, aes(x = harv, y = yield)) + geom_boxplot() + theme_bw(base_size = 12, base_family = "Arial")

a7 <- inner_join(a1.1, a6, by = "plot")
a7 <- a7 %>% spread(key = harv, value = yield)
head(a7)

a7$acc_num <- as.factor(a7$acc_num)
a7$rep <- as.factor(a7$rep)
a7$R <- as.factor(a7$row)
a7$C <- as.factor(a7$col)

nlevels(a7$acc_num)
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
