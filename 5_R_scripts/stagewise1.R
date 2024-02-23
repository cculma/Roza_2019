# stagewise
rm(list = ls())
library(StageWise)
library(SpATS)
library(tidyverse)

P2 <- read.csv("~/OneDrive - Washington State University (email.wsu.edu)/Roza_2019/yield/CSV_roza2019/Roza_2019_yield_20_21.csv")
head(P2)
colnames(P2)
# colnames(P2) <- c("cons","acc","id","col","row","rep","1","2","3","4","5","6","7","8","9","10")
colnames(P2)[3] <- "id"
P2 <- P2 %>% gather(key = "env", value = "yield", 7:16) %>% dplyr::select(7,3,6,4,5,8 )
head(P2)
dim(P2)
P2 <- na.omit(P2)
write.csv(P2, "~/OneDrive - Washington State University (email.wsu.edu)/Roza_2019/yield/CSV_roza2019/Roza_2019_Y.csv", quote = F, row.names = F)


directory <- "~/OneDrive - Washington State University (email.wsu.edu)/Roza_2019/yield/CSV_roza2019"
P5 <- file.path(directory, "Roza_2019_Y.csv")
P5
P3 <- read.csv("~/OneDrive - Washington State University (email.wsu.edu)/Roza_2019/yield/CSV_roza2019/Roza_2019_Y.csv")
head(P3)

effects <- data.frame(name=c("env", "col","row"),
                      fixed = c(F,F,F),
                      factor= c(F,T,T))
effects

model2 <- Stage1(filename=P5, 
                 traits="yield",
                 effects=effects, solver="spats", spline=c("col","row"))

stage1.vcov <- model2$vcov


predans <- asreml::predict.asreml(ans,classify="id",vcov = TRUE)
tmp <- predans$pvals[,c("id","predicted.value")]
colnames(tmp) <- c("id","BLUE")
vcov[[j]] <- predans$vcov
aic <- as.numeric(summary(ans)$aic)


summary(asreml_model)$bic
plot(model)
h2 <- vpredict(asreml_model, h2~v1 / (v1 + v2))

# obtaining BLUPs

BLUP <- summary(asreml_model, coef = T)$coef.random
