library(metan)
a1 <- inspect(data_ge2, verbose = FALSE)

setwd("~/Documents/git/Roza_2019/raw_data/")
a2 <- read.csv("Roza2019_yield_1.csv")
str(a2)
a2$acc_num <- as.factor(a2$acc_num)
a2$rep <- as.factor(a2$rep)
a3 <- inspect(a2, verbose = FALSE)

setwd("~/Documents/git/Roza_2019/pheno_data/")
write.csv(a3, "stats_metan.csv", quote = F, row.names = F)

