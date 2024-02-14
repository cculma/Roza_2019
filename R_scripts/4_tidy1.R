# tidy

library(tidyverse)

setwd("~/Documents/git/big_files/")
T1 <- read.csv("NS_DS.csv")
head(T1)
T1 <- T1 %>% gather(key = "env", value = "yield", 2:9) %>% separate(col = "env", into = c("Stress", "Year"), sep = "_", remove = T)
T1 <- T1 %>% relocate("yield", .before = "PC1")

write.csv(T1, "BLUE_Yi_sqrt_DArT_DS.csv", row.names = F, quote = F)
