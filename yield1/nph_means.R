rm(list = ls())
library(plyr)
library(tidyverse)
A1 <- read.csv("~/OneDrive - Washington State University (email.wsu.edu)/Roza_2019/yield/CSV_roza2019/Roza_2019_nph_20_21.csv")

str(A1)
col1 <- c("acc", "acc_num", "col", "row", "plant_num", "rep", "acc_rep")
A1[,col1] <- lapply(A1[,col1], factor)
A2 <- na.omit(A1)
A3 <- A2 %>% distinct(acc_rep, .keep_all = T)
A3 <- A3[,c(2:5,7,15)]
colnames(A3)
A3 <- A3 %>% remove_rownames() %>% column_to_rownames("acc_rep")

r1 <- ddply(A2, .(acc_rep), summarize, mean=mean(temp))
r2 <- ddply(A2, .(acc_rep), summarize, mean=mean(height_2019))
r3 <- ddply(A2, .(acc_rep), summarize, mean=mean(height_2020))
r4 <- cbind(r1, r2,r3)
r4 <- r4[,c(2,4,6)]
colnames(r4) <- c("temp", "height_2019", "height_2020")
rownames(r4) <- r1$acc_rep

A4 <- merge(as.data.frame(A3), as.data.frame(r4), by = "row.names", all = T)
write.csv(A4, "~/OneDrive - Washington State University (email.wsu.edu)/Roza_2019/yield/CSV_roza2019/Roza_2019_nph_temp.csv", quote = F, row.names = F)
