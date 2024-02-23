# spread pheno data
setwd("~/Documents/git/Roza_2019/pheno_data/")
s1 <- read.csv("spatial.csv")
s1 <- s1[,c(1:3)]
s1 <- s1 %>% spread(key = col, value = Plant_ID)

write.csv(s1, "spatial1.csv", row.names = F, quote = F)
