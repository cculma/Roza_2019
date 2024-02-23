rm(list = ls())
library(tidyverse)

setwd("~/Documents/git/Roza_2019/Roza_2019_Mr.Bean/")

a1 <- read.csv("Effects_mrbean.csv")

lev1 <- c("may_20","jun_20","jul_20","aug_20","sep_20",
          "may_21","jun_21","jul_21","aug_21","sep_21",
                   "jun_22","jul_22","aug_22","sep_22",
          "may_23","jun_23","jul_23","aug_23",
          "total_20","total_21","total_22")

lev2 <- c("may_20","jun_20","jul_20","aug_20","sep_20",
          "may_21","jun_21","jul_21","aug_21","sep_21",
          "jun_22","jul_22","aug_22","sep_22",
          "may_23","jun_23","jul_23","aug_23")


levels(a1$Trait)
a1$Trait <- as.factor(a1$Trait)
a1$Trait <- factor(a1$Trait, levels = lev1)
levels(a1$Trait)

a1.1 <- a1 %>% dplyr::select(1:3) %>% spread(key = Trait, value = predicted.values)
colnames(a1.1)

a1.2 <- a1 %>% dplyr::filter(!Trait %in% c("total_20","total_21","total_22")) 

ggplot(a1.2, aes(x=Trait, y=predicted.values)) + geom_boxplot()


a2 <- read.csv("predictions_2stage_ASReml_mrbean.csv")

levels(a2$trial)
a2$trial <- as.factor(a2$trial)
a2$trial <- factor(a2$trial, levels = lev2)
levels(a2$trial)

a2.1 <- a2 %>% dplyr::select(1:3) %>% spread(key = trial, value = predicted.value)

ggplot(a2, aes(x=trial, y=predicted.value)) + geom_boxplot()

a3 <- read.csv("Overall_predictions_2stage.csv")
a3.1 <- a3 %>% dplyr::select(1:2)
colnames(a3.1)[2] <- "all_yield"

colnames(a1.1)[1] <- "gen"
colnames(a1.1)
colnames(a2.1)

a3.2 <- left_join(a2.1, a3.1, by = "gen")

# write.csv(a1.1, "~/Documents/git/Roza_2019/Roza_2019_Mr.Bean/BLUP_ST1.csv", quote = F, row.names = F)
# write.csv(a3.2, "~/Documents/git/Roza_2019/Roza_2019_Mr.Bean/BLUP_ST2.csv", quote = F, row.names = F)

# pearson cor

a1.3 <- a1.1
colnames(a1.3)[2:22] <- paste0("ST1_", colnames(a1.3)[2:22])
colnames(a1.3)

colnames(a3.2)[2:20] <- paste0("ST2_", colnames(a3.2)[2:20])
colnames(a3.2)
a3.4 <- inner_join(a1.3, a3.2, by = "gen") %>% column_to_rownames("gen")

cor1 <- cor(a3.4, use = "complete.obs")

write.csv(cor1, "~/Documents/git/Roza_2019/Roza_2019_Mr.Bean/cor1.csv", quote = F, row.names = T)


