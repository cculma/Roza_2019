library(ggplot2)
library(tidyverse)
library(asreml)

# load("~/OneDrive - Washington State University (email.wsu.edu)/Roza_2019/GWASpoly_results/data_3_80177_year.RData")


# may 1
# jun 2
# jul 3
# aug 4
# sep 5
pheno <- read.csv("~/Documents/Cesar/blup_data/Roza2019/Analysis_2021/GWAS/pheno_nph.csv", row.names = 1)
head(pheno)
p1 <- pheno[1:(length(colnames(pheno))-3)]

p1 <- p1 %>% rownames_to_column(var = "gen")
colnames(p1)

p2 <- p1[,c(1:11)]

colnames(p2) <-  c("gen", "ST1_Yi_2020_1","ST1_Yi_2020_2", "ST1_Yi_2020_3", "ST1_Yi_2020_4", "ST1_Yi_2020_5", "ST1_Yi_2021_1", "ST1_Yi_2021_2", "ST1_Yi_2021_3", "ST1_Yi_2021_4", "ST1_Yi_2021_5")
p2 <- p2 %>% gather(key = "trait1", value = "BLUE", 2:11) %>% separate(2, into = c("stage", "trait", "year", "cut"), sep = "_", remove = F, convert = T, extra = "warn")
head(p2)
str(p2)
L2 <- colnames(p2)[1:(length(colnames(p2))-1)]
p2[,L2] <- lapply(p2[,L2], factor)

# ggplot(data = p2, aes(x = cut, y = BLUE, group = gen, color = gen)) + geom_line() + facet_grid(. ~ year) + theme(legend.position = "none")

trait_2 <- colnames(p2)[1:(length(colnames(p2))-1)]
p2[,trait_2] <- lapply(p2[,trait_2], factor)
str(p2)
hist(p2$BLUE)


########
# ASReml

# variance model for the residual term
model_uniform <- asreml(BLUE ~ year*cut,
                        residual = ~gen:cor(cut), data = p2)

summary(model_uniform)$bic

# the first autoregressive model (ar1) is the simplest model that allows correlations between pairs of measurements t decrease as the time between them increases.
# ar1 can only be fitted to equally spaced measurements.

# antedependence model
model_ante <- asreml(BLU3 ~ Year + cut + Year:cut,
                        residual = ~gen:ante(cut), data = p2)
summary(model_ante)$bic

wald(model_ante)

# predict tratment by time means
predict(model_ante, classify = "Year:cut")


P1 <- read.csv("~/Documents/Cesar/git/Norberg_2020/BLUE_values/split_data/ID_2018_1.csv")
head(P1)
P1 <- P1[,c(1:10,14)]
# P1 <- na.omit(P1)
model_uniform <- asreml::asreml(fixed = resp ~ 1 + gen, 
               random = ~+block + gen, residual = ~id(row):ar1(col), 						
               data = df, 
               na.action = list(x = "include", y = "include")) 	

model_uniform <- asreml::asreml(fixed = Yield ~1 + at(check, "control"):Treatment,
                                random = ~ + Block + at(check, "test"):Treatment, 
                                residual = ~id(row):id(col),
                                data = P1, na.action = list(x = "include", y = "include"))
str(P1)
check <- c(201, 202)

asreml::asreml(fixed = resp ~ 1 + gen, random = ~+block, residual = ~id(row):id(col), 
               data = df, na.action = list(x = "include", y = "include"))


model_uniform <- asreml(BLUE ~ year*cut,
                        residual = ~gen:cor(cut), data = p2)