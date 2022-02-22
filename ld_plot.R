# generate manhattan plot by trait

rm(list = ls())
# library(plotly)
# library(dplyr)
# library(ggpubr)
library(GWASpoly)
library(ggplot2)

# load("~/OneDrive - Washington State University (email.wsu.edu)/Roza_2019/GWASpoly_results/data_3_80177_year.RData")
load("data_3_80177_year.RData")
load("~/Documents/Cesar/blup_data/Roza2019/git/Roza2019/myplot1.RData")
myplot1 <- LD.plot(data_5) + theme_classic(base_family = "Arial", base_size = 12) + ggtitle("LD plot")
myplot2 <- myplot1 + theme_classic(base_family = "Arial", base_size = 12) + ggtitle("LD plot")


ggsave(filename = "~/Documents/Cesar/blup_data/Roza2019/git/Roza2019/myplot1.pdf", plot = myplot1, width = 6, height = 6)

ggsave(filename = "myplot1.jpg", plot = myplot1, width = 6, height = 6)
save(myplot1, file = "myplot1.RData")
save(myplot1, file = "myplot1.RData")