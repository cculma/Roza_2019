# tidy

library(tidyverse)

setwd("~/Documents/git/big_files/")
T1 <- read.csv("NS_DS.csv")
head(T1)
T1 <- T1 %>% gather(key = "env", value = "yield", 2:9) %>% separate(col = "env", into = c("Stress", "Year"), sep = "_", remove = T)
T1 <- T1 %>% relocate("yield", .before = "PC1")

write.csv(T1, "BLUE_Yi_sqrt_DArT_DS.csv", row.names = F, quote = F)

T2 <- T1 %>% dplyr::select(Roza2019_VCF_ID, yield, Year, Stress)

T2 <- T1 %>% spread(key = "Year", value = "yield")
head(T2)
T2 <- T2 %>% relocate("PC3", .after = "23") %>% relocate("PC2", .after = "23") %>% relocate("PC1", .after = "23") %>% relocate("Stress", .after = "PC3") 


setwd("~/Documents/git/big_files/")
write.csv(T2, "BLUE_Yi_sqrt_DArT_DS1.csv", row.names = F, quote = F)
