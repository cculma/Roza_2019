library(tidyverse)

setwd("~/Documents/git/Roza_2019/Roza_2019_MB/")

a1 <- read.csv("Effects_mrbean.csv")
a1.1 <- a1 %>% dplyr::select(1:3) %>% spread(key = Trait, value = predicted.values)
colnames(a1.1)

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

a1.2 <- a1 %>% dplyr::filter(!Trait %in% c("total_20","total_21","total_22")) 

ggplot(a1.2, aes(x=Trait, y=predicted.values)) + geom_boxplot()


a2 <- read.csv("predictions_2stage_ASReml_mrbean.csv")
a2.1 <- a2 %>% dplyr::select(1:3) %>% spread(key = trial, value = predicted.value)

levels(a2$trial)
a2$trial <- as.factor(a2$trial)
a2$trial <- factor(a2$trial, levels = lev2)

ggplot(a2, aes(x=trial, y=predicted.value)) + geom_boxplot()

