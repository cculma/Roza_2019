library(tidyverse)

setwd("~/medin297@umn.edu - Google Drive/My Drive/Roza2019/")

a1 <- read.csv("Roza_ID1.csv")
a2 <- read.csv("Roza_ID2.csv")
head(a1)
head(a2)

str(a1)
a1$Plant_ID <- as.character(a1$Plant_ID)
str(a2)

a3 <- right_join(a1, a2, by = "Plant_ID")
head(a3)
a4 <- a3 %>% distinct(Plant_ID, .keep_all = TRUE)

write.csv(a4, "Roza_ID3.csv", quote = F, row.names = F)
