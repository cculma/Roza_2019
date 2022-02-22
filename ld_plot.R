# generate manhattan plot by trait

rm(list = ls())
# library(plotly)
# library(dplyr)
# library(ggpubr)
library(GWASpoly)
library(ggplot2)

# load("~/OneDrive - Washington State University (email.wsu.edu)/Roza_2019/GWASpoly_results/data_3_80177_year.RData")
load("data_3_80177_year.RData")

myplot1 <- LD.plot(data_5) + theme_classic(base_family = "Arial", base_size = 12) + ggtitle("LD plot")

ggsave(filename = "myplot1.jpg", plot = myplot1, width = 6, height = 6)
save(myplot1, file = "myplot1.RData")