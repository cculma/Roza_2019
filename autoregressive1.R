library(ggplot2)
library(tidyverse)
load("~/OneDrive - Washington State University (email.wsu.edu)/Roza_2019/GWASpoly_results/data_3_80177_year.RData")
colnames(pheno)

# may 1
# jun 2
# jul 3
# aug 4
# sep 5

p1 <- pheno %>% dplyr::select(-c(24:26)) %>% rownames_to_column(var = "gen")
colnames(p1)

p2 <- p1[,c(1:11)]

colnames(p2) <-  c("gen", "ST1_Yi_2020_1","ST1_Yi_2020_2", "ST1_Yi_2020_3", "ST1_Yi_2020_4", "ST1_Yi_2020_5", "ST1_Yi_2021_1", "ST1_Yi_2021_2", "ST1_Yi_2021_3", "ST1_Yi_2021_4", "ST1_Yi_2021_5")
p2 <- p2 %>% gather(key = "trait1", value = "BLUE", 2:11) %>% separate(2, into = c("stage", "trait", "year", "cut"), sep = "_", remove = F, convert = T, extra = "warn")

ggplot(data = p2, aes(x = cut, y = BLUE, group = gen, color = gen)) + geom_line() + facet_grid(. ~ year) + theme(legend.position = "none")

trait_2 <- colnames(p2)[1:(length(colnames(p2))-1)]
p2[,trait_2] <- lapply(p2[,trait_2], factor)
str(p2)
hist(p2$BLUE)


########
# ASReml

# variance model for the residual term
model_uniform <- asreml(BLU3 ~ Year + cut + Year:cut,
                        residual = ~gen:cor(cut), data = p2)

summary(model_uniform)$bic

# the first autoregressive model (ar1) is the simplest model that allows correlations between pairs of measurements t decrease as the time between them increases.
# ar1 can only be fitted to equally spaced measurements.

# antedepence model
model_ante <- asreml(BLU3 ~ Year + cut + Year:cut,
                        residual = ~gen:ante(cut), data = p2)
summary(model_ante)$bic

wald(model_ante)

# predict tratment by time means
predict(model_ante, classify = "Year:cut")
